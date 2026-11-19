"""Keyword classes for gene products and Pfam descriptions (shared by survey.py and pfam_summary.py)."""
import re

# (class, case-insensitive words, case-sensitive gene symbols); first match wins
CLASSES = [
    ('CRISPR-Cas', r'CRISPR', r'\b(Cas\d+|Csm\d|Cmr\d|Csx\d*|Csy\d|Cse\d|Csn\d)'),
    ('phage/MGE', r'(?<!macro)phage|capsid|terminase|portal protein|\btail\b|baseplate|integrase|transposase|recombinase|'
                  r'excisionase|insertion element|conjuga|relaxase|plasmid|mobiliz|holin|endolysin|virion|'
                  r'antirepressor|reverse transcriptase', r'\b(IS\d+|Tra[A-Z]|Mob[A-Z]|Trb[A-Z])\b'),
    ('defence/TA', r'restriction|modification|DNA methyltransferase|DNA adenine methylase|toxin|abortive|'
                   r'Gabija|Zorya|Thoeris|BREX|DISARM|Septu|Druantia|Hachiman|PIN domain',
                   r'\b(Abi[A-Z]|HigA|HigB|RelE|RelB|ParE|VapC|VapB|MazF|HipA|Doc|Hsd[MRS])\b'),
    ('motility/chemotaxis/pili', r'flagell|chemotaxis|methyl-accepting|pilus|pilin|gliding|type IV pil|motility',
                                 r'\b(Flg[A-Z]|Fli[A-Z]|Flh[A-Z]|Che[A-Z]|Pil[A-Z]|Gld[A-Z]|Spr[A-Z]|Mot[AB])\b'),
    ('secretion', r'secretion|type VI|outer membrane protein assembly', r'\b(Tss[A-Z]|Hcp|VgrG|PorV|TolC|HlyD)\b'),
    ('cell surface/polysaccharide', r'glycosyltransferase|glycosyl transferase|O-antigen|lipopolysaccharide|'
                                    r'capsul|polysaccharide|epimerase|dehydratase|nucleotidyltransferase|'
                                    r'flippase|exopolysacchar|teichoic|sugar transferase|dTDP|UDP-|GDP-|CMP-',
                                    r'\b(Wza|Wzx|Wzy|Wzz|Wzc)\b'),
    ('CAZymes/glycan use', r'glycoside hydrolase|glycosyl hydrolase|sulfatase|TonB-dependent|alginate|lyase|'
                           r'glucanase|amylase|chitin|xylan|laminarin|fucosidase|galactosidase|glucosidase|'
                           r'mannosidase|arabinosidase|pectate|carbohydrate', r'\b(SusC|SusD)\b'),
    ('transport', r'transporter|permease|\bABC\b|porin|channel|efflux|symporter|antiporter|\bTRAP\b|'
                  r'uptake|exporter|pump', r'\bMFS\b'),
    ('housekeeping', r'ribosom|tRNA|release factor|DNA polymerase|RNA polymerase subunit|DNA-directed RNA polymerase|'
                     r'elongation factor|initiation factor|gyrase|topoisomerase|chaperon|DNA primase|'
                     r'replicative|DNA ligase|recombination protein Rec|DNA repair|ATP synthase', None),
    ('regulation/signalling', r'regulator|helix-turn-helix|\bHTH|sigma|histidine kinase|diguanylate|GGDEF|\bEAL\b|HD-GYP|transcription', None),
    ('metabolism/enzyme', r'(?<!rele)ase\b|dehydrogenase|reductase|oxidase|synth|kinase|transferase|hydrolase|oxygenase|'
                          r'isomerase|mutase|ligase|cytochrome', None),
    ('hypothetical/DUF', r'hypothetical|DUF|uncharacteri|domain-containing protein|unknown', None),
]


def fclass(prod):
    if not prod or prod == '?':
        return 'no hit'
    for name, ci, cs in CLASSES:
        if re.search(ci, prod, re.I) or (cs and re.search(cs, prod)):
            return name
    return 'other'
