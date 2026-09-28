#!/usr/bin/env python3
"""Export human WikiPathways interactions, map genes to Ensembl, and filter for brain expression.
Python 3.10+. Install: python -m pip install pandas openpyxl
See the accompanying README for filtering and output semantics.
"""
import argparse
import csv
import hashlib
import itertools
import json
import re
import sys
import time
import unicodedata
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlencode, unquote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

WP = 'http://vocabularies.wikipathways.org/wp#'
RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
RDFS = 'http://www.w3.org/2000/01/rdf-schema#'
DC = 'http://purl.org/dc/elements/1.1/'
DT = 'http://purl.org/dc/terms/'
PREFIX = f'PREFIX wp: <{WP}> PREFIX dcterms: <{DT}>\n'
COLUMNS = ['source_label', 'source_identifier', 'source_id_database',
           'target_label', 'target_identifier', 'target_id_database',
           'interaction_type', 'pathway_id', 'pathway_count', 'pathway_name',
           'evidence_count', 'source_database', 'organism', 'source_uri',
           'target_uri', 'interaction_uri', 'pathway_uri', 'directionality',
           'participant_uris', 'evidence_uris', 'source_count', 'target_count']

# Endpoints matched by the metabolite exclusion list (filled by organize(); used in checkpoints).
DELETED_NODES = set()

def joined(values):
    return ';'.join(sorted(set(values)))

def iri(value):
    if any(c in value for c in '<>"{}|^`\\\n\r '):
        raise ValueError(f'Unsafe IRI: {value!r}')
    return '<' + value + '>'

class Client:
    def __init__(self, endpoint, page_size=2000):
        self.endpoint, self.page_size = endpoint, page_size

    def query(self, query):
        url = self.endpoint + '?' + urlencode({'query': query, 'format': 'json'})
        for attempt in range(5):
            try:
                request = Request(url, headers={'Accept': 'application/sparql-results+json',
                                               'User-Agent': 'WikiPathwaysInteractionExport/1.0'})
                with urlopen(request, timeout=180) as response:
                    if response.headers.get('X-SQL-State', '00000') != '00000':
                        raise RuntimeError('Endpoint reported incomplete results')
                    data = json.load(response)
                return [{k: v['value'] for k, v in row.items()}
                        for row in data['results']['bindings']]
            except (HTTPError, URLError, TimeoutError, OSError, RuntimeError, ValueError) as exc:
                if isinstance(exc, HTTPError) and exc.code not in {408, 429, 500, 502, 503, 504}:
                    raise
                if attempt == 4:
                    raise
                print(f'Query failed ({type(exc).__name__}: {exc}); retry {attempt + 1}/4', file=sys.stderr)
                time.sleep(2 ** attempt)

    def pages(self, query, order):
        offset = 0
        size = self.page_size
        while True:
            try:
                page = self.query(f'{query}\nORDER BY {order}\nLIMIT {size} OFFSET {offset}')
            except (HTTPError, URLError, TimeoutError, OSError, RuntimeError) as exc:
                if isinstance(exc, HTTPError) and exc.code not in {408, 429, 500, 502, 503, 504}:
                    raise
                if size <= 100:
                    raise
                size = max(100, size // 2)
                print(f'Retrying offset {offset} with smaller page size {size}', file=sys.stderr)
                continue
            yield from page
            # Continue even after a short page: servers may impose a lower row cap.
            if not page:
                break
            offset += len(page)

def pathway_query(organism):
    return PREFIX + '''SELECT DISTINCT ?pathway WHERE {
      ?pathway a wp:Pathway ; wp:organismName ?organism .
      FILTER(STR(?organism) = ''' + json.dumps(organism) + ''')
    }'''

def triples_query(pathway):
    return PREFIX + '''SELECT DISTINCT ?s ?p ?o WHERE {
      { BIND(''' + iri(pathway) + ''' AS ?s) ?s ?p ?o }
      UNION {
        ?s a wp:Interaction ; dcterms:isPartOf ''' + iri(pathway) + ''' ; ?p ?o .
      }
      UNION {
        ?interaction a wp:Interaction ; dcterms:isPartOf ''' + iri(pathway) + ''' .
        VALUES ?role { wp:source wp:target wp:participants }
        ?interaction ?role ?s .
        ?s ?p ?o .
        FILTER(?p = <''' + RDFS + '''label> || STRSTARTS(STR(?p), CONCAT(STR(wp:), "bdb")))
      }
    }'''

def identifier_and_database(uri):
    """Use the endpoint's own RDF identifier, rather than a cross-mapping."""
    if not uri:
        return '', ''
    match = re.fullmatch(r'https?://identifiers\.org/([^/:]+)[/:](.+)', uri, re.I)
    if match:
        namespace, identifier = match.groups()
        return unquote(identifier), namespace.lower()
    if uri.startswith(('http://rdf.wikipathways.org/', 'https://rdf.wikipathways.org/')):
        # Keep the full identifier: a local graph ID alone is not globally unique.
        return uri, 'WikiPathways RDF'
    return uri, ''

def normalized_id(database, identifier):
    database = database.casefold()
    database = {'pubchem': 'pubchem.compound', 'chembl': 'chembl.compound'}.get(database, database)
    identifier = str(identifier).strip().upper()
    if database == 'hmdb' and re.fullmatch(r'HMDB\d+', identifier):
        identifier = 'HMDB' + str(int(identifier[4:]))
    if database == 'chembl.compound':
        identifier = identifier.replace('CHEMBL:', 'CHEMBL')
    return database, identifier

def normalized_label(value):
    return ' '.join(unicodedata.normalize('NFKC', str(value)).casefold().split())

class MetaboliteExclusions:
    def __init__(self, records):
        self.ids, self.labels = set(), set()
        for record in records:
            if record.get('chemicalLabel'):
                self.labels.add(normalized_label(record['chemicalLabel']))
            for column, database in [('wikiid', 'wikidata'), ('chembl_id', 'chembl.compound'),
                                     ('pubchem_id', 'pubchem.compound'), ('HMDB', 'hmdb')]:
                value = record.get(column)
                if value is not None and str(value).strip():
                    if isinstance(value, float) and value.is_integer():
                        value = int(value)
                    self.ids.add(normalized_id(database, value))

    @classmethod
    def load(cls, path):
        required = {'wikiid', 'chemicalLabel', 'chembl_id', 'pubchem_id', 'HMDB'}
        if path.suffix.lower() == '.xlsx':
            try:
                import openpyxl
            except ImportError as exc:
                raise RuntimeError('Reading XLSX requires: pip install openpyxl') from exc
            workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
            try:
                records = []
                found = False
                for sheet in workbook:
                    rows = iter(sheet.values)
                    headers = [str(v).strip() if v is not None else '' for v in next(rows, ())]
                    if required.issubset(headers):
                        found = True
                        records.extend(dict(zip(headers, row)) for row in rows)
                if not found:
                    raise ValueError('No worksheet has the expected metabolite-list headers')
            finally:
                workbook.close()
        else:
            with path.open(encoding='utf-8-sig', newline='') as handle:
                reader = csv.DictReader(handle)
                if not required.issubset(reader.fieldnames or []):
                    raise ValueError('Metabolite CSV is missing required headers')
                records = list(reader)
        result = cls(records)
        if not result.ids and not result.labels:
            raise ValueError('Metabolite exclusion list is empty')
        return result

    def matches(self, node, attrs):
        if not node:
            return False
        for label in attrs.get(RDFS + 'label', ()):
            if normalized_label(label) in self.labels:
                return True
        candidates = {node}
        for predicate, values in attrs.items():
            if predicate.startswith(WP + 'bdb'):
                candidates.update(values)
        for uri in candidates:
            identifier, database = identifier_and_database(uri)
            if normalized_id(database, identifier) in self.ids:
                return True
            match = re.fullmatch(r'https?://www.wikidata.org/entity/(Q\d+)', uri)
            if match and normalized_id('wikidata', match[1]) in self.ids:
                return True
        return False

def organize(pathway, triples, exclusions=None, removed=None):
    graph = defaultdict(lambda: defaultdict(set))
    for row in triples:
        graph[row['s']][row['p']].add(row['o'])
    meta = graph[pathway]
    identifiers = meta[DT + 'identifier'] | meta[DC + 'identifier'] | {pathway}
    ids = {m.group(0) for v in identifiers for m in re.finditer(r'\bWP\d+(?=_|\b)', v)}
    pathway_id = joined(ids) or pathway
    for interaction, attrs in sorted(list(graph.items())):
        if WP + 'Interaction' not in attrs[RDF + 'type']:
            continue
        sources, targets = attrs[WP + 'source'], attrs[WP + 'target']
        participants = attrs[WP + 'participants'] | sources | targets
        # Remove the whole RDF interaction so projections cannot retain a banned participant.
        matched = [node for node in participants if exclusions and exclusions.matches(node, graph[node])]
        if matched:
            DELETED_NODES.update(matched)
            if removed is not None:
                removed.add((pathway, interaction))
            continue
        types = {t.removeprefix(WP) for t in attrs[RDF + 'type'] if t.startswith(WP)}
        specific = types - {'Interaction', 'DirectedInteraction'}
        types = specific or (types - {'Interaction'}) or {'Interaction'}
        # Directly linked references only. Do not inherit pathway bibliography.
        evidence = attrs[DT + 'references'] | attrs[DC + 'references']
        if sources and targets:
            pairs = itertools.product(sorted(sources), sorted(targets))
            direction = 'directed'
        else:
            # Preserve participant-only and incomplete interactions without inventing edges.
            pairs = itertools.product(sorted(sources) or [''], sorted(targets) or [''])
            direction = 'incomplete' if sources or targets else 'undirected_or_unspecified'
        for source, target in pairs:
            def labels(node):
                return joined(graph[node][RDFS + 'label']) if node else ''
            source_id, source_db = identifier_and_database(source)
            target_id, target_db = identifier_and_database(target)
            yield dict(zip(COLUMNS, [
                labels(source), source_id, source_db, labels(target), target_id, target_db,
                joined(types), pathway_id, 0, joined(meta[DC + 'title'] | meta[DT + 'title']),
                len(evidence), 'WikiPathways', joined(meta[WP + 'organismName']), source,
                target, interaction, pathway, direction, joined(participants),
                joined(evidence), len(sources), len(targets)]))

def edge_key(row):
    # Use RDF identity, not potentially ambiguous labels or many-to-many mappings.
    endpoints = (row['source_uri'], row['target_uri'])
    if not all(endpoints):
        endpoints += (row['participant_uris'],)
    return (row['organism'], row['directionality'], row['interaction_type'], *endpoints)

def export_raw_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--organism', default='Homo sapiens')
    parser.add_argument('--endpoint', default='https://sparql.wikipathways.org/sparql')
    parser.add_argument('--output', type=Path, default=Path('human_wikipathways_interactions.csv'))
    parser.add_argument('--cache', type=Path, default=Path('wikipathways_cache'))
    parser.add_argument('--refresh', action='store_true', help='Replace cached pathway triples')
    parser.add_argument('--pathway', action='append', help='Optional WP identifier; repeatable')
    parser.add_argument('--page-size', type=int, default=2000)
    parser.add_argument('--exclude-metabolites', type=Path,
                        default=Path(__file__).with_name('metabolites to delete.xlsx'),
                        help='Metabolite list CSV or XLSX (default: metabolites to delete.xlsx beside the script)')
    args = parser.parse_args()
    if args.page_size < 1:
        parser.error('--page-size must be positive')
    exclusions = MetaboliteExclusions.load(args.exclude_metabolites)
    removed = set()
    client = Client(args.endpoint, args.page_size)
    pathways = sorted({r['pathway'] for r in client.pages(pathway_query(args.organism), '?pathway')})
    if args.pathway:
        requested = set(args.pathway)
        pathways = [p for p in pathways if requested.intersection(re.findall(r'\bWP\d+(?=_|\b)', p))]
    if not pathways:
        raise RuntimeError('No matching pathways returned; no output written.')
    args.cache.mkdir(parents=True, exist_ok=True)
    rows = []
    failures = []
    for index, pathway in enumerate(pathways, 1):
        print(f'[{index}/{len(pathways)}] {pathway}', file=sys.stderr)
        cache = args.cache / (hashlib.sha256(('all-bdb-v2' + args.endpoint + pathway).encode()).hexdigest() + '.json')
        if cache.exists() and not args.refresh:
            triples = json.loads(cache.read_text(encoding='utf-8'))
        else:
            try:
                triples = list(client.pages(triples_query(pathway), '?s ?p ?o'))
            except (HTTPError, URLError, TimeoutError, OSError, RuntimeError, ValueError) as exc:
                failures.append({'pathway_uri': pathway, 'error': f'{type(exc).__name__}: {exc}'})
                print(f'FAILED {pathway}: {exc}; continuing with other pathways', file=sys.stderr)
                continue
            temporary = cache.with_suffix('.tmp')
            temporary.write_text(json.dumps(triples, ensure_ascii=False), encoding='utf-8')
            temporary.replace(cache)
        rows.extend(organize(pathway, triples, exclusions, removed))
    rows = [row for row in rows if any(str(row.get(field) or '').strip() for field in
            ('source_label', 'source_identifier', 'source_id_database', 'target_label', 'target_identifier'))]
    excluded_databases = {'wikipathways', 'wikipathways rdf', 'aop.events', 'go'}
    rows = [row for row in rows if not any(
        str(row.get(field) or '').strip().casefold().lstrip(':') in excluded_databases
        for field in ('source_id_database', 'target_id_database'))]
    memberships = defaultdict(set)
    for row in rows:
        memberships[edge_key(row)].add(row['pathway_id'])
    for row in rows:
        row['pathway_count'] = len(memberships[edge_key(row)])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    destination = args.output.with_name(args.output.stem + '.partial' + args.output.suffix) if failures else args.output
    temporary = destination.with_suffix(destination.suffix + '.tmp')
    with temporary.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(destination)
    report = args.output.with_name(args.output.stem + '.status.json')
    report_tmp = report.with_suffix('.json.tmp')
    report_tmp.write_text(json.dumps({
        'complete': not failures, 'output': str(destination),
        'pathways_requested': len(pathways), 'pathways_completed': len(pathways) - len(failures),
        'failed_pathways': failures, 'rows': len(rows),
        'note': 'Pathway counts cover successfully retrieved pathways only.'
    }, indent=2), encoding='utf-8')
    report_tmp.replace(report)
    print(f'Wrote {len(rows):,} rows from {len(pathways) - len(failures):,} pathway resources to {destination}')
    print(f'Excluded {len(removed):,} pathway/interaction records using {args.exclude_metabolites}')
    if failures:
        print(f'INCOMPLETE: {len(failures)} failed pathways listed in {report}. '
              'Run the same command again without --refresh to retry missing pathways.', file=sys.stderr)
        return 2
    return 0



import pandas as pd
from datetime import datetime, timezone

RDF_TYPE = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type'

CATEGORIES = ['Protein–protein', 'Metabolite–protein', 'Metabolite–metabolite', 'Other/unknown']

EXCLUDED = {'wikipathways', 'wikipathways rdf', 'aop.events'}

CHEMICAL_DBS = {'chebi', 'hmdb', 'pubchem.compound', 'chembl.compound', 'kegg.compound',
                'lipidmaps', 'cas', 'chemspider', 'drugbank', 'pubchem.substance'}

GENE_DBS = {'ncbigene', 'ensembl', 'uniprot', 'hgnc', 'hgnc.symbol', 'refseq', 'ncbiprotein',
            'kegg.genes'}

GENE_ID_SCOPES = {
    'ncbigene': 'entrezgene',
    'hgnc': 'hgnc',
    'hgnc.symbol': 'symbol',
    'uniprot': 'uniprot',
    'refseq': 'refseq',
    'ncbiprotein': 'refseq',
}

KEGG_GENE_DBS = {'kegg.genes'}

MYGENE_URL = 'https://mygene.info/v3/query'

MYGENE_BATCH = 1000

HPA_BRAIN_URL = 'https://www.proteinatlas.org/download/tsv/rna_brain_region_hpa.tsv.zip'

def load_brain_expression(path, threshold=1.0, refresh=False, download=False):
    """Full HPA human regional nTPM table; select genes with max regional nTPM >= threshold."""
    if not 0 < threshold < float('inf'):
        raise ValueError('Brain threshold must be a positive finite nTPM value')
    if download and (refresh or not path.exists()):
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + '.download')
        for attempt in range(5):
            try:
                request = Request(HPA_BRAIN_URL, headers={'User-Agent': 'BrainInteractionIntegration/1.0'})
                with urlopen(request, timeout=180) as response, temporary.open('wb') as target:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        target.write(chunk)
                # Check schema before accepting a download into the cache.
                header = pd.read_csv(temporary, sep='\t', compression='zip', nrows=0)
                if not {'Gene', 'Brain region', 'nTPM'}.issubset(header.columns):
                    raise ValueError('Downloaded HPA data lacks Gene/Brain region/nTPM columns')
                temporary.replace(path)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
    if not path.exists():
        raise FileNotFoundError(f'Brain-expression file not found: {path}')
    expression = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, compression='infer')
    required = {'Gene', 'Brain region', 'nTPM'}
    if not required.issubset(expression.columns):
        raise ValueError('Brain expression input must be a regional TSV/TSV.ZIP with Gene, Brain region, nTPM columns')
    expression['ensembl_id'] = expression.Gene.map(gene_id)
    if expression.ensembl_id.eq('').any() or expression['Brain region'].str.strip().eq('').any():
        raise ValueError('Brain expression data contains invalid gene IDs or empty regions')
    values = pd.to_numeric(expression.nTPM, errors='coerce')
    invalid = (expression.nTPM.str.strip().ne('') & values.isna()) | values.lt(0) | values.eq(float('inf'))
    if invalid.any():
        raise ValueError('Brain expression data contains invalid nTPM values')
    expression['nTPM'] = values
    # One value per gene/region; repeated records do not inflate region counts.
    regional = expression.groupby(['ensembl_id', 'Brain region'], as_index=False).nTPM.max()
    all_genes = regional.groupby('ensembl_id').agg(max_brain_nTPM=('nTPM', 'max'),
                                                  measured_regions=('nTPM', 'count'))
    passed = regional[regional.nTPM.ge(threshold)].groupby('ensembl_id')['Brain region'].agg(
        lambda v: ';'.join(sorted(set(v))))
    all_genes['expressed_regions'] = passed.reindex(all_genes.index).fillna('')
    all_genes['brain_expressed'] = all_genes.max_brain_nTPM.ge(threshold)
    if 'Gene name' in expression:
        symbols = expression.groupby('ensembl_id')['Gene name'].agg(lambda v: ';'.join(sorted(set(v) - {''})))
        all_genes['gene_symbol'] = symbols.reindex(all_genes.index).fillna('')
    selected = set(all_genes.index[all_genes.brain_expressed])
    if not selected:
        raise ValueError('No genes pass the selected brain-expression threshold')
    details = {'source_url': HPA_BRAIN_URL if download else None, 'input_file': str(path.resolve()),
               'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'processed_at_utc': datetime.now(timezone.utc).isoformat(),
               'threshold_nTPM': threshold, 'rule': 'nTPM >= threshold in at least one region',
               'regions': sorted(set(regional['Brain region'])), 'genes_in_table': len(all_genes),
               'brain_expressed_genes': len(selected),
               'note': 'Brain-expressed, not brain-exclusive; includes all regions supplied by HPA, including spinal cord.'}
    return selected, all_genes.reset_index(), details

def split(value):
    return sorted({x.strip() for x in str(value or '').split(';') if x.strip()})

def gene_id(value):
    value = unquote(str(value)).rsplit('/', 1)[-1].removeprefix('ensembl:').upper()
    return value.split('.')[0] if re.fullmatch(r'ENSG\d+(?:\.\d+)?', value) else ''

def identity(uri, database='', identifier=''):
    match = re.fullmatch(r'https?://identifiers\.org/([^/:]+)[/:](.+)', uri)
    if match:
        database, identifier = match[1], unquote(match[2])
    database = database.strip().casefold().lstrip(':')
    database = {'pubchem': 'pubchem.compound', 'chembl': 'chembl.compound'}.get(database, database)
    identifier = identifier.strip()
    if database == 'hgnc':
        identifier = re.sub(r'^HGNC:', '', identifier, flags=re.I)
    if database == 'hmdb' and re.fullmatch(r'HMDB\d+', identifier, re.I):
        identifier = 'HMDB' + str(int(identifier[4:])).zfill(7)
    if database == 'chembl.compound':
        identifier = identifier.upper().replace('CHEMBL:', 'CHEMBL')
    if database == 'ensembl' and gene_id(identifier):
        identifier = gene_id(identifier)
    return database, identifier

def endpoint_uri(row, side):
    uri = row.get(side + '_uri', '').strip()
    if uri:
        return uri
    db, identifier = identity('', row.get(side + '_id_database', ''), row.get(side + '_identifier', ''))
    if db and identifier and ';' not in identifier:
        return f'https://identifiers.org/{db}/{identifier}'
    return ''

def sparql(query, endpoint):
    for attempt in range(5):
        try:
            req = Request(endpoint + '?' + urlencode({'query': query, 'format': 'json'}),
                          headers={'Accept': 'application/sparql-results+json',
                                   'User-Agent': 'HumanInteractionIntegration/1.0'})
            with urlopen(req, timeout=180) as response:
                if response.headers.get('X-SQL-State', '00000') != '00000':
                    raise RuntimeError('SPARQL returned an incomplete result')
                return json.load(response)['results']['bindings']
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)

def mapping_query(nodes):
    for node in nodes:
        if any(c in node for c in '<>"{}|^`\\\r\n '):
            raise ValueError(f'Invalid endpoint URI: {node!r}')
    values = ' '.join('<' + n + '>' for n in nodes)
    return PREFIX + f'''SELECT DISTINCT ?node ?predicate ?value WHERE {{
      VALUES ?node {{ {values} }}
      ?node dcterms:isPartOf ?pathway .
      ?pathway a wp:Pathway ; wp:organismName ?species .
      FILTER(STR(?species) = "Homo sapiens")
      VALUES ?predicate {{ <{RDF_TYPE}> wp:bdbEnsembl }}
      ?node ?predicate ?value .
    }}'''

def fetch_annotations(nodes, cache, endpoint, refresh=False):
    """Cache only complete batches; failed queries abort rather than become 'unmapped'."""
    cache.mkdir(parents=True, exist_ok=True)
    result = {node: {'genes': set(), 'types': set()} for node in nodes}
    nodes = sorted(nodes)
    for start in range(0, len(nodes), 40):
        batch = nodes[start:start+40]
        query = mapping_query(batch)
        file = cache / (hashlib.sha256((endpoint + query).encode()).hexdigest() + '.json')
        if file.exists() and not refresh:
            bindings = json.loads(file.read_text(encoding='utf-8'))
        else:
            bindings, offset = [], 0
            while True:
                page = sparql(query + f' ORDER BY ?node ?predicate ?value LIMIT 500 OFFSET {offset}', endpoint)
                if not page:
                    break
                bindings.extend(page)
                offset += len(page)
            temporary = file.with_suffix('.tmp')
            temporary.write_text(json.dumps(bindings), encoding='utf-8')
            temporary.replace(file)
        for binding in bindings:
            node, pred, value = (binding[k]['value'] for k in ('node', 'predicate', 'value'))
            if pred == RDF_TYPE:
                result[node]['types'].add(value)
            elif gene_id(value):
                result[node]['genes'].add(gene_id(value))
        if start % 400 == 0 or start + 40 >= len(nodes):
            print(f'Mapped/typed {min(start+40, len(nodes)):,}/{len(nodes):,} WikiPathways endpoints', flush=True)
    return result

def classify_wp_node(db, identifier, annotation):
    """Return (is_chemical, is_gene, genes) for a WikiPathways endpoint, from its RDF
    annotation alone (before any external gene-ID mapping is consulted)."""
    types = {t.rsplit('#', 1)[-1] for t in annotation['types']}
    genes = set(annotation['genes'])
    if db == 'ensembl' and gene_id(identifier):
        genes.add(gene_id(identifier))
    is_chemical = 'Metabolite' in types or db in CHEMICAL_DBS
    is_gene = bool(types & {'GeneProduct', 'Protein', 'Dna', 'Rna', 'Gene'}) or db in GENE_DBS or bool(genes)
    return is_chemical, is_gene, genes

def kegg_query_term(identifier):
    """KEGG human gene IDs are 'hsa:<NCBI Gene ID>'; the numeric part is an Entrez Gene ID."""
    term = identifier.rsplit(':', 1)[-1].strip()
    return term if term.isdigit() else ''

def pending_gene_lookups(wp_df, annotations):
    """(db, identifier) pairs for gene-type WikiPathways endpoints that have no Ensembl
    xref in the WikiPathways RDF and therefore need an external ID-mapping lookup."""
    pending = set()
    for row in wp_df.to_dict('records'):
        for side in ('source', 'target'):
            uri = endpoint_uri(row, side)
            if not uri:
                continue
            db, identifier = identity(uri, row.get(side + '_id_database', ''), row.get(side + '_identifier', ''))
            if db in EXCLUDED or uri.startswith(('http://rdf.wikipathways.org/', 'https://rdf.wikipathways.org/')):
                continue
            annotation = annotations.get(uri, {'genes': set(), 'types': set()})
            is_chemical, is_gene, genes = classify_wp_node(db, identifier, annotation)
            if is_gene and not is_chemical and not genes and db and identifier:
                pending.add((db, identifier))
    return pending

def mygene_query(terms_by_scope, cache, refresh=False, species='human'):
    """Batch-query mygene.info's /v3/query for one or more scopes.
    terms_by_scope: {scope: [query terms]}. Returns {(scope, query_term): {'ensembl': set(), 'symbol': set()}}.
    Only complete batches are cached, matching fetch_annotations' caching contract."""
    cache.mkdir(parents=True, exist_ok=True)
    result = {}
    for scope, raw_terms in terms_by_scope.items():
        terms = sorted({t for t in raw_terms if t})
        for start in range(0, len(terms), MYGENE_BATCH):
            batch = terms[start:start + MYGENE_BATCH]
            key = hashlib.sha256((MYGENE_URL + '|' + species + '|' + scope + '|' + '|'.join(batch)).encode()).hexdigest()
            file = cache / f'{key}.json'
            if file.exists() and not refresh:
                hits = json.loads(file.read_text(encoding='utf-8'))
            else:
                payload = urlencode({'q': ','.join(batch), 'scopes': scope,
                                      'fields': 'symbol,ensembl.gene,taxid', 'species': species}).encode()
                request = Request(MYGENE_URL, data=payload,
                                   headers={'Content-Type': 'application/x-www-form-urlencoded',
                                            'User-Agent': 'HumanInteractionIntegration/1.0'})
                hits = None
                for attempt in range(5):
                    try:
                        with urlopen(request, timeout=120) as response:
                            hits = json.load(response)
                        break
                    except Exception:
                        if attempt == 4:
                            raise
                        time.sleep(2 ** attempt)
                temporary = file.with_suffix('.tmp')
                temporary.write_text(json.dumps(hits), encoding='utf-8')
                temporary.replace(file)
                time.sleep(0.34)  # be polite to the public API (~3 requests/second)
            for hit in hits:
                if hit.get('notfound') or hit.get('taxid', 9606) != 9606:
                    continue
                entry = result.setdefault((scope, hit.get('query', '')), {'ensembl': set(), 'symbol': set()})
                ensembl_field = hit.get('ensembl')
                items = ensembl_field if isinstance(ensembl_field, list) else [ensembl_field] if ensembl_field else []
                for item in items:
                    gid = gene_id((item or {}).get('gene', ''))
                    if gid:
                        entry['ensembl'].add(gid)
                if hit.get('symbol'):
                    entry['symbol'].add(hit['symbol'])
            print(f'mygene.info [{scope}]: queried {min(start+MYGENE_BATCH, len(terms)):,}/{len(terms):,} terms', flush=True)
    return result

def map_wp_genes_to_ensembl(pending, cache, refresh=False):
    """pending: set of (db, identifier) WikiPathways gene endpoints without a direct Ensembl
    xref. Returns {(db, identifier): {'ensembl': set(), 'symbol': set()}} for everything that
    could be resolved via mygene.info. Pairs whose db has no configured scope are skipped
    (left to fall back to the old gene_product_without_ensembl behaviour) and reported."""
    terms_by_scope = defaultdict(list)
    term_lookup = defaultdict(list)
    skipped = defaultdict(int)
    for db, identifier in pending:
        if db in KEGG_GENE_DBS:
            scope, term = 'entrezgene', kegg_query_term(identifier)
        else:
            scope, term = GENE_ID_SCOPES.get(db), identifier
        if not scope or not term:
            skipped[db] += 1
            continue
        terms_by_scope[scope].append(term)
        term_lookup[(scope, term)].append((db, identifier))
    if skipped:
        detail = ', '.join(f'{db} ({n:,})' for db, n in sorted(skipped.items()))
        print(f'No external ID-mapping scope configured for: {detail}; these stay as their own db:identifier node.', flush=True)
    if not terms_by_scope:
        return {}
    hits = mygene_query(terms_by_scope, cache, refresh)
    result = {}
    for (scope, term), entry in hits.items():
        for db, identifier in term_lookup.get((scope, term), []):
            result[(db, identifier)] = entry
    return result

def ensembl_symbols(ensembl_ids, cache, refresh=False):
    """Bulk-resolve canonical gene symbols for a set of Ensembl gene IDs via mygene.info."""
    if not ensembl_ids:
        return {}
    hits = mygene_query({'ensembl.gene': sorted(ensembl_ids)}, cache, refresh)
    symbols = {}
    for (scope, query), entry in hits.items():
        if entry['symbol']:
            symbols[query] = sorted(entry['symbol'])[0]
    return symbols

def resolve_endpoint(row, side, annotations, external_map):
    uri = endpoint_uri(row, side)
    db, identifier = identity(uri, row.get(side+'_id_database', ''), row.get(side+'_identifier', ''))
    if db in EXCLUDED or uri.startswith(('http://rdf.wikipathways.org/', 'https://rdf.wikipathways.org/')):
        return [], 'other', 'excluded_database'
    annotation = annotations.get(uri, {'genes': set(), 'types': set()})
    is_chemical, is_gene, genes = classify_wp_node(db, identifier, annotation)
    if is_chemical and is_gene:
        kind, status = 'other', 'conflicting_type_annotations'
    elif is_chemical:
        kind, status = 'metabolite', 'metabolite_identifier_retained'
    elif is_gene:
        status_tag = 'mapped'
        if not genes:
            mapped = external_map.get((db, identifier))
            if mapped and mapped['ensembl']:
                genes = set(mapped['ensembl'])
                status_tag = 'mapped_via_external_id_mapping'
        if genes:
            return sorted(genes), 'protein', status_tag if len(genes) == 1 else status_tag + '_multiple'
        kind, status = 'protein', 'gene_product_without_ensembl'
    else:
        kind, status = 'other', 'unknown_type'
    key = f'{db}:{identifier}' if db and identifier else uri
    return ([key] if key else []), kind, status if key else 'missing_endpoint'

def interaction_category(a, b):
    if a == b == 'protein':
        return CATEGORIES[0]
    if {a, b} == {'protein', 'metabolite'}:
        return CATEGORIES[1]
    if a == b == 'metabolite':
        return CATEGORIES[2]
    return CATEGORIES[3]

def save_csv(path, rows, fields):
    """Always write headers, including empty audit tables; replace output atomically."""
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def brain_status(node, kind, selected, measured):
    if kind == 'metabolite':
        return 'metabolite_exempt'
    if kind != 'protein':
        return 'unclassified_endpoint'
    gid = gene_id(node)
    if not gid:
        return 'unmapped_gene'
    if gid in selected:
        return 'brain_expressed'
    return 'below_threshold' if gid in measured else 'no_expression_data'


def print_checkpoints(table, audit, retained):
    """Report unique endpoint -> Ensembl -> brain counts; no filtering changes."""
    directed = table.directionality.eq('directed')
    print('\n===== INTERACTIONS -> ENSEMBL -> BRAIN EXPRESSION =====')
    print(f'[1] Retrieved records after the existing raw-export filters: {len(table):,}')
    print(f'    Unique interaction URIs: {table.interaction_uri.nunique():,}')
    print(f'    Directed records: {int(directed.sum()):,}; '
          f'incomplete: {int(table.directionality.eq("incomplete").sum()):,}; '
          f'participant-only/unspecified: '
          f'{int(table.directionality.eq("undirected_or_unspecified").sum()):,}')

    scopes = [('ALL source/target endpoints', audit),
              ('DIRECTED records only', [a for a in audit
                                        if a.get('directionality') == 'directed'])]
    for scope, rows in scopes:
        products, metabolites, others = set(), set(), set()
        mapped_products, brain_products = set(), set()
        ensembl, brain_ensembl, below, no_data = set(), set(), set(), set()
        for a in rows:
            original = a.get('original_uri', '')
            if not original:
                continue
            kind = a.get('node_type', 'other')
            if kind == 'metabolite':
                metabolites.add(a.get('node_id') or original)
                continue
            if kind != 'protein':
                others.add(original)
                continue
            products.add(original)
            gid = gene_id(a.get('node_id', ''))
            if not gid:
                continue
            mapped_products.add(original)
            ensembl.add(gid)
            status = a.get('brain_filter_status', '')
            if status == 'brain_expressed':
                brain_products.add(original)
                brain_ensembl.add(gid)
            elif status == 'below_threshold':
                below.add(gid)
            elif status == 'no_expression_data':
                no_data.add(gid)

        print(f'\n    Scope: {scope}')
        print(f'[2] BEFORE Ensembl mapping: {len(products):,} unique gene/protein endpoints; '
              f'{len(metabolites):,} unique metabolite identifiers; '
              f'{len(others):,} other/unclassified endpoints.')
        print(f'[3] ENSEMBL MAPPING: {len(mapped_products):,} of {len(products):,} '
              f'gene/protein endpoints mapped; {len(products - mapped_products):,} unmapped.')
        print(f'    Mapping yields {len(ensembl):,} UNIQUE Ensembl gene IDs.')
        print(f'[4] BRAIN FILTER: {len(brain_ensembl):,} of {len(ensembl):,} '
              f'unique Ensembl gene IDs pass the configured brain-expression threshold.')
        print(f'    Below threshold: {len(below - brain_ensembl):,}; '
              f'no expression data: {len(no_data - brain_ensembl - below):,}; '
              f'other/unreported status: {len(ensembl - brain_ensembl - below - no_data):,}.')
        print(f'    These passing genes represent {len(brain_products):,} original '
              f'gene/protein endpoints with at least one passing mapping.')

    final_genes = {gene_id(r[s + '_id']) for r in retained for s in ('source', 'target')
                   if r[s + '_type'] == 'protein'} - {''}
    final_metabolites = {r[s + '_id'] for r in retained for s in ('source', 'target')
                        if r[s + '_type'] == 'metabolite'}
    print(f'\n[5] FINAL retained interaction network: {len(final_genes):,} unique Ensembl genes '
          f'and {len(final_metabolites):,} unique metabolite identifiers.')
    print('    Final-network counts can be lower: a passing gene also needs a retained interaction partner.')
    print(f'    Exclusion-list matches: {len(DELETED_NODES):,} unique endpoints.')
    print('Counting: source/target endpoints only; repeated rows are deduplicated. '
          'Unassigned participant-only nodes are not included in endpoint counts.')
    print('Before mapping, gene/protein endpoints are identified by original URI, not unique biological gene. '
          'Multiple identifiers can map to one gene, and one endpoint can map to several genes.')
    print('=====================================================\n')



def run(args):
    output = args.output_folder
    output.mkdir(parents=True, exist_ok=True)
    status_file = output / 'export.status.json'
    status_file.write_text(json.dumps({'complete': False, 'state': 'running'}), encoding='utf-8')
    brain_path = args.brain_expression or args.cache / 'rna_brain_region_hpa.tsv.zip'
    brain_genes, brain_summary, brain_details = load_brain_expression(
        brain_path, args.brain_threshold, args.refresh_brain_expression,
        args.brain_expression is None)
    measured = set(brain_summary.loc[brain_summary.measured_regions.gt(0), 'ensembl_id'])
    brain_summary.to_csv(output / 'brain_expression_gene_summary.csv', index=False)
    (output / 'brain_expression_provenance.json').write_text(json.dumps(brain_details, indent=2), encoding='utf-8')
    raw_path = output / 'wikipathways_raw_interactions.csv'
    command = ['export_wp_brain.py', '--organism', 'Homo sapiens',
               '--endpoint', args.endpoint, '--output', str(raw_path),
               '--cache', str(args.cache / 'pathways'), '--page-size', str(args.page_size),
               '--exclude-metabolites', str(args.exclude_metabolites)]
    if args.refresh:
        command.append('--refresh')
    for pathway in args.pathway or []:
        command.extend(['--pathway', pathway])
    old_argv = sys.argv
    try:
        sys.argv = command
        code = export_raw_main()
    finally:
        sys.argv = old_argv
    if code:
        raise RuntimeError('WikiPathways retrieval was incomplete. See the raw export status and rerun without --refresh.')
    table = pd.read_csv(raw_path, dtype=str, keep_default_na=False)
    if not table.organism.str.strip().eq('Homo sapiens').all():
        raise ValueError('Non-human or unspecified organism found in raw export')
    records = table.to_dict('records')
    uris = {endpoint_uri(r, side) for r in records for side in ('source', 'target')} - {''}
    annotations = fetch_annotations(uris, args.cache / 'annotations', args.endpoint, args.refresh)
    pending = pending_gene_lookups(table, annotations)
    external = map_wp_genes_to_ensembl(pending, args.cache / 'mygene', args.refresh)
    retained, rejected, audit = [], [], []
    ok = {'brain_expressed', 'metabolite_exempt'}
    for number, row in enumerate(records, 2):
        ends, kinds = [], []
        for side in ('source', 'target'):
            ids, kind, mapping = resolve_endpoint(row, side, annotations, external)
            ends.append(ids)
            kinds.append(kind)
            for node in ids or ['']:
                audit.append({'input_csv_row': number, 'side': side,
                              'directionality': row['directionality'],
                              'original_uri': endpoint_uri(row, side),
                              'original_label': row.get(side + '_label', ''),
                              'node_id': node, 'node_type': kind, 'mapping_status': mapping,
                              'brain_filter_status': brain_status(node, kind, brain_genes, measured)})
        if not all(ends):
            rejected.append({'input_csv_row': number, 'source_id': ';'.join(ends[0]),
                             'target_id': ';'.join(ends[1]), 'reason': 'missing_or_excluded_endpoint',
                             'interaction_uri': row['interaction_uri']})
            continue
        for a, b in itertools.product(*ends):
            sa, sb = (brain_status(n, k, brain_genes, measured) for n, k in zip((a, b), kinds))
            if sa not in ok or sb not in ok or a == b:
                rejected.append({'input_csv_row': number, 'source_id': a, 'target_id': b,
                                 'reason': 'self_loop' if a == b else f'source={sa};target={sb}',
                                 'interaction_uri': row['interaction_uri']})
                continue
            retained.append({**row, 'source_original_label': row['source_label'],
                             'target_original_label': row['target_label'],
                             'source_id': a, 'target_id': b,
                             'source_type': kinds[0], 'target_type': kinds[1],
                             'source_ensembl': gene_id(a), 'target_ensembl': gene_id(b),
                             'source_brain_status': sa, 'target_brain_status': sb,
                             'interaction_category': interaction_category(*kinds)})
    audit_fields = ['input_csv_row', 'side', 'original_uri', 'original_label', 'node_id',
                    'node_type', 'mapping_status', 'brain_filter_status']
    save_csv(output / 'endpoint_audit.csv', audit, audit_fields)
    save_csv(output / 'removed_interactions.csv', rejected,
             ['input_csv_row', 'source_id', 'target_id', 'reason', 'interaction_uri'])
    print_checkpoints(table, audit, retained)
    genes = {r[s + '_ensembl'] for r in retained for s in ('source', 'target')} - {''}
    symbols = ensembl_symbols(genes, args.cache / 'mygene', args.refresh)
    for r in retained:
        for side in ('source', 'target'):
            # Never reuse a potentially multi-gene source label as a canonical gene symbol.
            if r[side + '_type'] == 'protein':
                r[side + '_label'] = symbols.get(r[side + '_id'], r[side + '_id'])
    extra_fields = ['source_id', 'source_type', 'source_ensembl', 'source_original_label',
                    'target_id', 'target_type', 'target_ensembl', 'target_original_label',
                    'source_brain_status', 'target_brain_status', 'interaction_category']
    save_csv(output / 'wp_brain_interactions.csv', retained, COLUMNS + extra_fields)
    # Additional explicitly undirected, unique-pair table for network tools.
    nodes, unique = {}, {}
    for r in retained:
        for side in ('source', 'target'):
            node = r[side + '_id']
            item = nodes.setdefault(node, {'labels': set(), 'type': r[side + '_type']})
            if item['type'] != r[side + '_type']:
                raise ValueError(f'Conflicting endpoint types for {node}')
            item['labels'].add(r[side + '_label'])
        key = tuple(sorted((r['source_id'], r['target_id'])))
        item = unique.setdefault(key, {'pathways': set(), 'types': set()})
        item['pathways'].update(split(r['pathway_id']))
        item['types'].update(split(r['interaction_type']))
    edge_rows = []
    for (a, b), item in sorted(unique.items()):
        edge_rows.append({'source_id': a, 'source_label': joined(nodes[a]['labels']),
                         'source_type': nodes[a]['type'], 'source_ensembl': gene_id(a),
                         'target_id': b, 'target_label': joined(nodes[b]['labels']),
                         'target_type': nodes[b]['type'], 'target_ensembl': gene_id(b),
                         'interaction_category': interaction_category(nodes[a]['type'], nodes[b]['type']),
                         'database_count': 1, 'source_database': 'WikiPathways',
                         'pathway_id': joined(item['pathways']), 'pathway_count': len(item['pathways']),
                         'interaction_types': joined(item['types'])})
    edge_fields = ['source_id', 'source_label', 'source_type', 'source_ensembl',
                   'target_id', 'target_label', 'target_type', 'target_ensembl',
                   'interaction_category', 'database_count', 'source_database',
                   'pathway_id', 'pathway_count', 'interaction_types']
    save_csv(output / 'wp_brain_edges.csv', edge_rows, edge_fields)
    status_file.write_text(json.dumps({'complete': True, 'organism': 'Homo sapiens',
        'raw_records': len(records), 'retained_interaction_records': len(retained),
        'unique_undirected_edges': len(edge_rows), 'removed_pair_records': len(rejected),
        'brain_threshold_nTPM': args.brain_threshold, 'brain_expressed_genes': len(brain_genes),
        'metabolites_exempt': True, 'finished_utc': datetime.now(timezone.utc).isoformat(),
        'note': 'wp_brain_interactions.csv preserves direction; wp_brain_edges.csv is undirected. '
                'Multiple gene mappings are expanded and individually filtered; inspect endpoint_audit.csv.'},
        indent=2), encoding='utf-8')
    print(f'Finished: {len(retained):,} interaction records; {len(edge_rows):,} unique undirected edges. {output}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-folder', type=Path, default=Path('wp_brain_output'))
    parser.add_argument('--cache', type=Path, default=Path('wp_brain_cache'))
    parser.add_argument('--endpoint', default='https://sparql.wikipathways.org/sparql')
    parser.add_argument('--pathway', action='append', help='Optional WP ID; repeat for several pathways')
    parser.add_argument('--page-size', type=int, default=2000)
    parser.add_argument('--refresh', action='store_true', help='Refresh pathway and identifier caches')
    parser.add_argument('--exclude-metabolites', type=Path,
                        default=Path(__file__).with_name('metabolites to delete.xlsx'))
    parser.add_argument('--brain-expression', type=Path,
                        help='Local HPA regional TSV/TSV.ZIP; otherwise download the HPA table')
    parser.add_argument('--brain-threshold', type=float, default=1.0)
    parser.add_argument('--refresh-brain-expression', action='store_true')
    args = parser.parse_args()
    if args.page_size < 1:
        parser.error('--page-size must be positive')
    try:
        run(args)
    except Exception as exc:
        args.output_folder.mkdir(parents=True, exist_ok=True)
        (args.output_folder / 'export.status.json').write_text(json.dumps({
            'complete': False, 'error': f'{type(exc).__name__}: {exc}',
            'note': 'Do not use outputs from this folder until a run reports complete=true.'}, indent=2), encoding='utf-8')
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
