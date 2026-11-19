"""Pick PanGBank GTDB_all 1.0.0 (GTDB r226) pangenomes for the data-size experiments.

  A (level 2, context model vs number of genomes): species stratified by genome count
    (15-29, 30-59, 60-119, 120-249, 250-499, 500-999, 1000-MAXG), round-robin over phyla inside each stratum
    so that no phylum dominates; N_PER_BIN species per stratum.
  B (level 1, genus shrinkage): genera with >= 3 PanGBank species, round-robin over phyla; up to
    MAX_SP_PER_GENUS species per genus (smallest genome counts first, so downloads stay small).

Writes into OUTDIR:
  catalog.tsv          all GTDB_all 1.0.0 pangenomes with taxonomy and size (cached; delete to refresh)
  species_A.tsv        selection for A
  genera_B.tsv         selection for B (one row per species, grouped by genus)
  tasks.tsv            one line per worker task: 'species<TAB>id' or 'genus<TAB>genus_name<TAB>id,id,...',
                       sorted by expected cost, most expensive first
  downloads.tsv        id, name, genomes, md5 of every .h5 needed

usage: select_species.py OUTDIR [N_PER_BIN=45] [N_GENERA=40] [MAXG=1500] [SEED=1]
"""
import json, os, sys, time, random, urllib.parse, urllib.request
from collections import defaultdict

API = "https://pangbank-api.genoscope.cns.fr"
out = sys.argv[1]
N_PER_BIN = int(sys.argv[2]) if len(sys.argv) > 2 else 45
N_GENERA = int(sys.argv[3]) if len(sys.argv) > 3 else 40
MAXG = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 1
MAX_SP_PER_GENUS = 8
MAXG_B = 300             # genus members above this many genomes are skipped (download size)
BINS = [(15, 29), (30, 59), (60, 119), (120, 249), (250, 499), (500, 999), (1000, MAXG)]
os.makedirs(out, exist_ok=True)


def get(path, **q):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(f"{API}{path}?{urllib.parse.urlencode(q)}", timeout=120) as r:
                return json.load(r)
        except Exception as e:
            print(f"retry {attempt}: {e}", file=sys.stderr); time.sleep(10 * (attempt + 1))
    raise RuntimeError(path)


cat_f = f"{out}/catalog.tsv"
cols = ["id", "name", "genomes", "families", "persistent", "shell", "cloud", "modules", "mean_compl", "md5",
        "domain", "phylum", "class", "order", "family", "genus", "species"]
if not os.path.exists(cat_f):
    rows, off = [], 0
    while True:
        page = get("/pangenomes/", collection_name="GTDB_all", release_version="1.0.0",
                   only_latest_release="false", limit=100, offset=off)
        for h in page:
            tx = {t["rank"]: t["name"][3:] for t in h["taxonomy"]["taxa"]}
            rows.append([h["id"], h["name"], h["genome_count"], h["family_count"], h["persistent_family_count"],
                         h["shell_family_count"], h["cloud_family_count"], h["module_count"],
                         round(h["mean_completeness"], 1), h["file_md5sum"]] +
                        [tx.get(r, "") for r in ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
        off += len(page)
        if len(page) < 100: break
        time.sleep(0.3)
    with open(cat_f, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")
cat = [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in open(cat_f).readlines()[1:]]
for c in cat: c["genomes"] = int(c["genomes"])
print(f"catalog: {len(cat)} pangenomes, {len({c['phylum'] for c in cat})} phyla", file=sys.stderr)
rng = random.Random(SEED)


def round_robin(items, key, n):
    """pick n items, cycling over groups (phyla) in random order"""
    groups = defaultdict(list)
    for it in items: groups[key(it)].append(it)
    for g in groups.values(): rng.shuffle(g)
    order = list(groups); rng.shuffle(order)
    picked = []
    while len(picked) < n and any(groups.values()):
        for g in order:
            if groups[g] and len(picked) < n: picked.append(groups[g].pop())
    return picked


# ---- A
selA = []
for lo, hi in BINS:
    pool = [c for c in cat if lo <= c["genomes"] <= hi]
    selA += round_robin(pool, lambda c: c["phylum"], N_PER_BIN)
selA.sort(key=lambda c: -c["genomes"])

# ---- B
by_genus = defaultdict(list)
for c in cat: by_genus[(c["phylum"], c["genus"])].append(c)
gen = [(k, v) for k, v in by_genus.items() if len(v) >= 3]
pickedG = round_robin(gen, lambda kv: kv[0][0], N_GENERA)
selB = []
for (phy, g), members in pickedG:
    members = sorted([m for m in members if m["genomes"] <= MAXG_B], key=lambda m: m["genomes"])[:MAX_SP_PER_GENUS]
    if len(members) >= 3: selB.append((g, members))

with open(f"{out}/species_A.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for c in selA: f.write("\t".join(str(c[k]) for k in cols) + "\n")
with open(f"{out}/genera_B.tsv", "w") as f:
    f.write("genus_group\t" + "\t".join(cols) + "\n")
    for g, ms in selB:
        for c in ms: f.write(g + "\t" + "\t".join(str(c[k]) for k in cols) + "\n")

# cost proxy: A ~ min(n, 800)^1.5; B ~ total genomes
tasks = [(min(c["genomes"], 800) ** 1.5 * 2, f"species\t{c['id']}") for c in selA]
tasks += [(sum(m["genomes"] for m in ms) * 20, f"genus\t{g}\t{','.join(str(m['id']) for m in ms)}") for g, ms in selB]
tasks.sort(key=lambda t: -t[0])
with open(f"{out}/tasks.tsv", "w") as f:
    for _, t in tasks: f.write(t + "\n")
need = {c["id"]: c for c in selA}
for _, ms in selB:
    for c in ms: need[c["id"]] = c
with open(f"{out}/downloads.tsv", "w") as f:
    for c in sorted(need.values(), key=lambda c: c["genomes"]):
        f.write(f"{c['id']}\t{c['name']}\t{c['genomes']}\t{c['md5']}\n")
tot = sum(c["genomes"] for c in need.values())
print(f"A: {len(selA)} species ({len({c['phylum'] for c in selA})} phyla); B: {len(selB)} genera, "
      f"{sum(len(m) for _, m in selB)} species ({len({m[0]['phylum'] for _, m in selB})} phyla); "
      f"{len(need)} .h5 files, {tot} genomes, ~{tot * 0.5 / 1000:.0f} GB at 0.5 MB/genome", file=sys.stderr)
