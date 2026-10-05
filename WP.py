#!/usr/bin/env python3
"""Export human WikiPathways interactions, map genes to Ensembl, keep brain-expressed genes.
Python 3.10+. Requires pandas (and openpyxl for an XLSX exclusion list)."""
import argparse
import csv
import hashlib
import itertools
import json
import re
import shutil
import sys
import time
import unicodedata
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import unquote, urlencode
from urllib.request import Request, urlopen

import pandas as pd

WP = 'http://vocabularies.wikipathways.org/wp#'
RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
RDFS = 'http://www.w3.org/2000/01/rdf-schema#'
DC = 'http://purl.org/dc/elements/1.1/'
DT = 'http://purl.org/dc/terms/'
RDF_TYPE = RDF + 'type'
PREFIX = f'PREFIX wp: <{WP}> PREFIX dcterms: <{DT}>\n'
WP_RDF_PREFIXES = ('http://rdf.wikipathways.org/', 'https://rdf.wikipathways.org/')
UNSAFE_IRI = '<>"{}|^`\\\n\r '
RETRY_CODES = {408, 429, 500, 502, 503, 504}

COLUMNS = ['source_label', 'source_identifier', 'source_id_database',
           'target_label', 'target_identifier', 'target_id_database',
           'interaction_type', 'pathway_id', 'pathway_count', 'pathway_name',
           'evidence_count', 'source_database', 'organism', 'source_uri',
           'target_uri', 'interaction_uri', 'pathway_uri', 'directionality',
           'participant_uris', 'evidence_uris', 'source_count', 'target_count']
CATEGORIES = ['Protein–protein', 'Metabolite–protein', 'Metabolite–metabolite', 'Other/unknown']

EXCLUDED = {'wikipathways', 'wikipathways rdf', 'aop.events'}
RAW_EXCLUDED = EXCLUDED | {'go'}
CHEMICAL_DBS = {'chebi', 'hmdb', 'pubchem.compound', 'chembl.compound', 'kegg.compound',
                'lipidmaps', 'cas', 'chemspider', 'drugbank', 'pubchem.substance'}
GENE_DBS = {'ncbigene', 'ensembl', 'uniprot', 'hgnc', 'hgnc.symbol', 'refseq', 'ncbiprotein', 'kegg.genes'}
# WikiPathways identifier database -> mygene.info query scope
GENE_ID_SCOPES = {'ncbigene': 'entrezgene', 'hgnc': 'hgnc', 'hgnc.symbol': 'symbol',
                  'uniprot': 'uniprot', 'refseq': 'refseq', 'ncbiprotein': 'refseq'}
MYGENE_URL = 'https://mygene.info/v3/query'
MYGENE_BATCH = 1000
HPA_BRAIN_URL = 'https://www.proteinatlas.org/download/tsv/rna_brain_region_hpa.tsv.zip'

DELETED_NODES = set()  # endpoints matched by the metabolite exclusion list


# ---------------------------------------------------------------- utilities

def joined(values):
    return ';'.join(sorted(set(values)))


def split(value):
    return sorted({x.strip() for x in str(value or '').split(';') if x.strip()})


def iri(value):
    if any(c in value for c in UNSAFE_IRI):
        raise ValueError(f'Unsafe IRI: {value!r}')
    return f'<{value}>'


def retry(fn, tries=5, http_only=False):
    """Call fn with exponential backoff. With http_only, fail fast on non-transient HTTP codes."""
    for attempt in range(tries):
        try:
            return fn()
        except (OSError, RuntimeError, ValueError) as exc:
            if (http_only and isinstance(exc, HTTPError) and exc.code not in RETRY_CODES) or attempt == tries - 1:
                raise
            print(f'Retry {attempt + 1}/{tries - 1} ({type(exc).__name__}: {exc})', file=sys.stderr)
            time.sleep(2 ** attempt)


def write_json(path, data, **kwargs):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, **kwargs), encoding='utf-8')
    tmp.replace(path)


def save_csv(path, rows, fields):
    """Atomic write; headers are always written, even for empty tables."""
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


class Client:
    def __init__(self, endpoint, page_size=2000):
        self.endpoint, self.page_size = endpoint, page_size

    def query(self, query):
        url = self.endpoint + '?' + urlencode({'query': query, 'format': 'json'})

        def fetch():
            request = Request(url, headers={'Accept': 'application/sparql-results+json',
                                            'User-Agent': 'WikiPathwaysInteractionExport/1.0'})
            with urlopen(request, timeout=180) as response:
                if response.headers.get('X-SQL-State', '00000') != '00000':
                    raise RuntimeError('Endpoint reported incomplete results')
                data = json.load(response)
            return [{k: v['value'] for k, v in row.items()} for row in data['results']['bindings']]

        return retry(fetch, http_only=True)

    def pages(self, query, order):
        offset, size = 0, self.page_size
        while True:
            try:
                page = self.query(f'{query}\nORDER BY {order}\nLIMIT {size} OFFSET {offset}')
            except (OSError, RuntimeError, ValueError) as exc:
                if (isinstance(exc, HTTPError) and exc.code not in RETRY_CODES) or size <= 100:
                    raise
                size = max(100, size // 2)
                print(f'Retrying offset {offset} with page size {size}', file=sys.stderr)
                continue
            if not page:  # a short page is not the end: servers may cap rows
                return
            yield from page
            offset += len(page)


# ---------------------------------------------------------------- identifiers

def identifier_and_database(uri):
    """Identifier and database as given by the endpoint's own RDF."""
    if not uri:
        return '', ''
    match = re.fullmatch(r'https?://identifiers\.org/([^/:]+)[/:](.+)', uri, re.I)
    if match:
        return unquote(match[2]), match[1].lower()
    if uri.startswith(WP_RDF_PREFIXES):
        return uri, 'WikiPathways RDF'  # local graph IDs are not globally unique
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
    return f'https://identifiers.org/{db}/{identifier}' if db and identifier and ';' not in identifier else ''


# ---------------------------------------------------------------- metabolite exclusions

class MetaboliteExclusions:
    COLUMNS = [('wikiid', 'wikidata'), ('chembl_id', 'chembl.compound'),
               ('pubchem_id', 'pubchem.compound'), ('HMDB', 'hmdb')]

    def __init__(self, records):
        self.ids, self.labels = set(), set()
        for record in records:
            if record.get('chemicalLabel'):
                self.labels.add(normalized_label(record['chemicalLabel']))
            for column, database in self.COLUMNS:
                value = record.get(column)
                if value is not None and str(value).strip():
                    if isinstance(value, float) and value.is_integer():
                        value = int(value)
                    self.ids.add(normalized_id(database, value))

    @classmethod
    def load(cls, path):
        required = {'wikiid', 'chemicalLabel', 'chembl_id', 'pubchem_id', 'HMDB'}
        if path.suffix.lower() == '.xlsx':
            import openpyxl
            workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
            records = []
            try:
                for sheet in workbook:
                    rows = iter(sheet.values)
                    headers = [str(v).strip() if v is not None else '' for v in next(rows, ())]
                    if required.issubset(headers):
                        records.extend(dict(zip(headers, row)) for row in rows)
            finally:
                workbook.close()
            if not records:
                raise ValueError('No worksheet has the expected metabolite-list headers')
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
        if any(normalized_label(label) in self.labels for label in attrs.get(RDFS + 'label', ())):
            return True
        candidates = {node}
        for predicate, values in attrs.items():
            if predicate.startswith(WP + 'bdb'):
                candidates.update(values)
        for uri in candidates:
            if normalized_id(*reversed(identifier_and_database(uri))) in self.ids:
                return True
            wikidata = re.fullmatch(r'https?://www.wikidata.org/entity/(Q\d+)', uri)
            if wikidata and normalized_id('wikidata', wikidata[1]) in self.ids:
                return True
        return False


# ---------------------------------------------------------------- WikiPathways retrieval

def pathway_query(organism):
    return PREFIX + ('SELECT DISTINCT ?pathway WHERE { ?pathway a wp:Pathway ; wp:organismName ?organism . '
                     f'FILTER(STR(?organism) = {json.dumps(organism)}) }}')


def triples_query(pathway):
    p = iri(pathway)
    return PREFIX + f'''SELECT DISTINCT ?s ?p ?o WHERE {{
      {{ BIND({p} AS ?s) ?s ?p ?o }}
      UNION {{ ?s a wp:Interaction ; dcterms:isPartOf {p} ; ?p ?o . }}
      UNION {{
        ?interaction a wp:Interaction ; dcterms:isPartOf {p} .
        VALUES ?role {{ wp:source wp:target wp:participants }}
        ?interaction ?role ?s .
        ?s ?p ?o .
        FILTER(?p = <{RDFS}label> || STRSTARTS(STR(?p), CONCAT(STR(wp:), "bdb")))
      }}
    }}'''


def organize(pathway, triples, exclusions=None, removed=None):
    graph = defaultdict(lambda: defaultdict(set))
    for row in triples:
        graph[row['s']][row['p']].add(row['o'])
    meta = graph[pathway]
    identifiers = meta[DT + 'identifier'] | meta[DC + 'identifier'] | {pathway}
    ids = {m[0] for v in identifiers for m in re.finditer(r'\bWP\d+(?=_|\b)', v)}
    pathway_id = joined(ids) or pathway
    pathway_name = joined(meta[DC + 'title'] | meta[DT + 'title'])
    organism = joined(meta[WP + 'organismName'])
    label = lambda node: joined(graph[node][RDFS + 'label']) if node else ''

    for interaction in sorted(graph):
        attrs = graph[interaction]
        if WP + 'Interaction' not in attrs[RDF_TYPE]:
            continue
        sources, targets = attrs[WP + 'source'], attrs[WP + 'target']
        participants = attrs[WP + 'participants'] | sources | targets
        # Drop the whole interaction so no projection can keep a banned participant.
        matched = [n for n in participants if exclusions and exclusions.matches(n, graph[n])]
        if matched:
            DELETED_NODES.update(matched)
            if removed is not None:
                removed.add((pathway, interaction))
            continue
        types = {t.removeprefix(WP) for t in attrs[RDF_TYPE] if t.startswith(WP)}
        types = (types - {'Interaction', 'DirectedInteraction'}) or (types - {'Interaction'}) or {'Interaction'}
        evidence = attrs[DT + 'references'] | attrs[DC + 'references']  # direct links only
        if sources and targets:
            direction = 'directed'
        else:
            direction = 'incomplete' if sources or targets else 'undirected_or_unspecified'
        for source, target in itertools.product(sorted(sources) or [''], sorted(targets) or ['']):
            source_id, source_db = identifier_and_database(source)
            target_id, target_db = identifier_and_database(target)
            yield dict(zip(COLUMNS, [
                label(source), source_id, source_db, label(target), target_id, target_db,
                joined(types), pathway_id, 0, pathway_name, len(evidence), 'WikiPathways',
                organism, source, target, interaction, pathway, direction,
                joined(participants), joined(evidence), len(sources), len(targets)]))


def edge_key(row):
    endpoints = (row['source_uri'], row['target_uri'])
    if not all(endpoints):
        endpoints += (row['participant_uris'],)
    return (row['organism'], row['directionality'], row['interaction_type'], *endpoints)


def export_raw(endpoint, output, cache, page_size, refresh, requested, exclusions_path, organism='Homo sapiens'):
    exclusions = MetaboliteExclusions.load(exclusions_path)
    removed, rows, failures = set(), [], []
    client = Client(endpoint, page_size)
    pathways = sorted({r['pathway'] for r in client.pages(pathway_query(organism), '?pathway')})
    if requested:
        pathways = [p for p in pathways if set(requested) & set(re.findall(r'\bWP\d+(?=_|\b)', p))]
    if not pathways:
        raise RuntimeError('No matching pathways returned; no output written.')
    cache.mkdir(parents=True, exist_ok=True)

    for index, pathway in enumerate(pathways, 1):
        print(f'[{index}/{len(pathways)}] {pathway}', file=sys.stderr)
        file = cache / (hashlib.sha256(('all-bdb-v2' + endpoint + pathway).encode()).hexdigest() + '.json')
        if file.exists() and not refresh:
            triples = json.loads(file.read_text(encoding='utf-8'))
        else:
            try:
                triples = list(client.pages(triples_query(pathway), '?s ?p ?o'))
            except (OSError, RuntimeError, ValueError) as exc:
                failures.append({'pathway_uri': pathway, 'error': f'{type(exc).__name__}: {exc}'})
                print(f'FAILED {pathway}: {exc}; continuing', file=sys.stderr)
                continue
            write_json(file, triples, ensure_ascii=False)
        rows.extend(organize(pathway, triples, exclusions, removed))

    has_endpoint = ('source_label', 'source_identifier', 'source_id_database', 'target_label', 'target_identifier')
    rows = [r for r in rows if any(str(r.get(f) or '').strip() for f in has_endpoint)]
    rows = [r for r in rows if not any(
        str(r.get(f) or '').strip().casefold().lstrip(':') in RAW_EXCLUDED
        for f in ('source_id_database', 'target_id_database'))]
    memberships = defaultdict(set)
    for r in rows:
        memberships[edge_key(r)].add(r['pathway_id'])
    for r in rows:
        r['pathway_count'] = len(memberships[edge_key(r)])

    output.parent.mkdir(parents=True, exist_ok=True)
    destination = output.with_name(output.stem + '.partial' + output.suffix) if failures else output
    save_csv(destination, rows, COLUMNS)
    report = output.with_name(output.stem + '.status.json')
    write_json(report, {
        'complete': not failures, 'output': str(destination),
        'pathways_requested': len(pathways), 'pathways_completed': len(pathways) - len(failures),
        'failed_pathways': failures, 'rows': len(rows),
        'note': 'Pathway counts cover successfully retrieved pathways only.'}, indent=2)
    print(f'Wrote {len(rows):,} rows from {len(pathways) - len(failures):,} pathways to {destination}')
    print(f'Excluded {len(removed):,} pathway/interaction records using {exclusions_path}')
    if failures:
        print(f'INCOMPLETE: {len(failures)} failed pathways listed in {report}. '
              'Rerun without --refresh to retry the missing ones.', file=sys.stderr)
        return 2
    return 0


# ---------------------------------------------------------------- brain expression

def load_brain_expression(path, threshold=1.0, refresh=False, download=False):
    """HPA regional table; a gene passes if max regional nTPM >= threshold."""
    if not 0 < threshold < float('inf'):
        raise ValueError('Brain threshold must be a positive finite nTPM value')
    need = {'Gene', 'Brain region', 'nTPM'}

    def fetch():
        tmp = path.with_suffix(path.suffix + '.download')
        request = Request(HPA_BRAIN_URL, headers={'User-Agent': 'BrainInteractionIntegration/1.0'})
        with urlopen(request, timeout=180) as response, tmp.open('wb') as target:
            shutil.copyfileobj(response, target)
        if not need.issubset(pd.read_csv(tmp, sep='\t', compression='zip', nrows=0).columns):
            raise ValueError('Downloaded HPA data lacks Gene/Brain region/nTPM columns')
        tmp.replace(path)

    if download and (refresh or not path.exists()):
        path.parent.mkdir(parents=True, exist_ok=True)
        retry(fetch)
    if not path.exists():
        raise FileNotFoundError(f'Brain-expression file not found: {path}')

    expr = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, compression='infer')
    if not need.issubset(expr.columns):
        raise ValueError('Brain expression input must be a regional TSV/TSV.ZIP with Gene, Brain region, nTPM columns')
    expr['ensembl_id'] = expr.Gene.map(gene_id)
    if expr.ensembl_id.eq('').any() or expr['Brain region'].str.strip().eq('').any():
        raise ValueError('Brain expression data contains invalid gene IDs or empty regions')
    values = pd.to_numeric(expr.nTPM, errors='coerce')
    if ((expr.nTPM.str.strip().ne('') & values.isna()) | values.lt(0) | values.eq(float('inf'))).any():
        raise ValueError('Brain expression data contains invalid nTPM values')
    expr['nTPM'] = values

    regional = expr.groupby(['ensembl_id', 'Brain region'], as_index=False).nTPM.max()
    genes = regional.groupby('ensembl_id').agg(max_brain_nTPM=('nTPM', 'max'),
                                              measured_regions=('nTPM', 'count'))
    passed = regional[regional.nTPM.ge(threshold)].groupby('ensembl_id')['Brain region'].agg(
        lambda v: ';'.join(sorted(set(v))))
    genes['expressed_regions'] = passed.reindex(genes.index).fillna('')
    genes['brain_expressed'] = genes.max_brain_nTPM.ge(threshold)
    if 'Gene name' in expr:
        names = expr.groupby('ensembl_id')['Gene name'].agg(lambda v: ';'.join(sorted(set(v) - {''})))
        genes['gene_symbol'] = names.reindex(genes.index).fillna('')
    selected = set(genes.index[genes.brain_expressed])
    if not selected:
        raise ValueError('No genes pass the selected brain-expression threshold')
    details = {'source_url': HPA_BRAIN_URL if download else None, 'input_file': str(path.resolve()),
               'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'processed_at_utc': datetime.now(timezone.utc).isoformat(),
               'threshold_nTPM': threshold, 'rule': 'nTPM >= threshold in at least one region',
               'regions': sorted(set(regional['Brain region'])), 'genes_in_table': len(genes),
               'brain_expressed_genes': len(selected),
               'note': 'Brain-expressed, not brain-exclusive; includes all HPA regions, spinal cord included.'}
    return selected, genes.reset_index(), details


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


# ---------------------------------------------------------------- Ensembl mapping

def mapping_query(nodes):
    values = ' '.join(iri(n) for n in nodes)
    return PREFIX + f'''SELECT DISTINCT ?node ?predicate ?value WHERE {{
      VALUES ?node {{ {values} }}
      ?node dcterms:isPartOf ?pathway .
      ?pathway a wp:Pathway ; wp:organismName ?species .
      FILTER(STR(?species) = "Homo sapiens")
      VALUES ?predicate {{ <{RDF_TYPE}> wp:bdbEnsembl }}
      ?node ?predicate ?value .
    }}'''


def fetch_annotations(nodes, cache, endpoint, refresh=False):
    """RDF types and Ensembl xrefs per endpoint. Only complete batches are cached."""
    cache.mkdir(parents=True, exist_ok=True)
    client = Client(endpoint, 500)
    result = {n: {'genes': set(), 'types': set()} for n in nodes}
    nodes = sorted(nodes)
    for start in range(0, len(nodes), 40):
        query = mapping_query(nodes[start:start + 40])
        file = cache / (hashlib.sha256(('v2' + endpoint + query).encode()).hexdigest() + '.json')
        if file.exists() and not refresh:
            rows = json.loads(file.read_text(encoding='utf-8'))
        else:
            rows = list(client.pages(query, '?node ?predicate ?value'))
            write_json(file, rows)
        for row in rows:
            if row['predicate'] == RDF_TYPE:
                result[row['node']]['types'].add(row['value'])
            elif gene_id(row['value']):
                result[row['node']]['genes'].add(gene_id(row['value']))
        if start % 400 == 0 or start + 40 >= len(nodes):
            print(f'Mapped/typed {min(start + 40, len(nodes)):,}/{len(nodes):,} WikiPathways endpoints', flush=True)
    return result


def classify_wp_node(db, identifier, annotation):
    """(is_chemical, is_gene, ensembl genes) from the RDF annotation alone."""
    types = {t.rsplit('#', 1)[-1] for t in annotation['types']}
    genes = set(annotation['genes'])
    if db == 'ensembl' and gene_id(identifier):
        genes.add(gene_id(identifier))
    is_chemical = 'Metabolite' in types or db in CHEMICAL_DBS
    is_gene = bool(types & {'GeneProduct', 'Protein', 'Dna', 'Rna', 'Gene'}) or db in GENE_DBS or bool(genes)
    return is_chemical, is_gene, genes


def pending_gene_lookups(table, annotations):
    """(db, identifier) of gene endpoints with no Ensembl xref in the RDF."""
    pending = set()
    for row in table.to_dict('records'):
        for side in ('source', 'target'):
            uri = endpoint_uri(row, side)
            if not uri or uri.startswith(WP_RDF_PREFIXES):
                continue
            db, identifier = identity(uri, row.get(side + '_id_database', ''), row.get(side + '_identifier', ''))
            if db in EXCLUDED:
                continue
            is_chemical, is_gene, genes = classify_wp_node(
                db, identifier, annotations.get(uri, {'genes': set(), 'types': set()}))
            if is_gene and not is_chemical and not genes and db and identifier:
                pending.add((db, identifier))
    return pending


def mygene_query(terms_by_scope, cache, refresh=False, species='human'):
    """Batched mygene.info lookups -> {(scope, term): {'ensembl': set, 'symbol': set}}."""
    cache.mkdir(parents=True, exist_ok=True)
    result = {}
    for scope, raw in terms_by_scope.items():
        terms = sorted({t for t in raw if t})
        for start in range(0, len(terms), MYGENE_BATCH):
            batch = terms[start:start + MYGENE_BATCH]
            key = hashlib.sha256('|'.join([MYGENE_URL, species, scope, *batch]).encode()).hexdigest()
            file = cache / f'{key}.json'
            if file.exists() and not refresh:
                hits = json.loads(file.read_text(encoding='utf-8'))
            else:
                payload = urlencode({'q': ','.join(batch), 'scopes': scope,
                                     'fields': 'symbol,ensembl.gene,taxid', 'species': species}).encode()
                request = Request(MYGENE_URL, data=payload,
                                  headers={'Content-Type': 'application/x-www-form-urlencoded',
                                           'User-Agent': 'HumanInteractionIntegration/1.0'})

                def post():
                    with urlopen(request, timeout=120) as response:
                        return json.load(response)

                hits = retry(post)
                write_json(file, hits)
                time.sleep(0.34)  # public API: stay near 3 requests/second
            for hit in hits:
                if hit.get('notfound') or hit.get('taxid', 9606) != 9606:
                    continue
                entry = result.setdefault((scope, hit.get('query', '')), {'ensembl': set(), 'symbol': set()})
                field = hit.get('ensembl')
                for item in field if isinstance(field, list) else [field] if field else []:
                    if gene_id((item or {}).get('gene', '')):
                        entry['ensembl'].add(gene_id(item['gene']))
                if hit.get('symbol'):
                    entry['symbol'].add(hit['symbol'])
            print(f'mygene.info [{scope}]: queried {min(start + MYGENE_BATCH, len(terms)):,}/{len(terms):,} terms', flush=True)
    return result


def map_wp_genes_to_ensembl(pending, cache, refresh=False):
    """{(db, identifier): {'ensembl', 'symbol'}} for pending gene endpoints that mygene.info resolves.
    Databases without a configured scope are reported and keep their own db:identifier node."""
    terms_by_scope, lookup, skipped = defaultdict(list), defaultdict(list), defaultdict(int)
    for db, identifier in pending:
        if db == 'kegg.genes':  # 'hsa:<Entrez ID>'
            scope, term = 'entrezgene', identifier.rsplit(':', 1)[-1].strip()
            term = term if term.isdigit() else ''
        else:
            scope, term = GENE_ID_SCOPES.get(db), identifier
        if not scope or not term:
            skipped[db] += 1
            continue
        terms_by_scope[scope].append(term)
        lookup[(scope, term)].append((db, identifier))
    if skipped:
        print('No ID-mapping scope for: ' + ', '.join(f'{d} ({n:,})' for d, n in sorted(skipped.items()))
              + '; these stay as db:identifier nodes.', flush=True)
    if not terms_by_scope:
        return {}
    return {pair: entry for key, entry in mygene_query(terms_by_scope, cache, refresh).items()
            for pair in lookup.get(key, [])}


def ensembl_symbols(ensembl_ids, cache, refresh=False):
    if not ensembl_ids:
        return {}
    hits = mygene_query({'ensembl.gene': sorted(ensembl_ids)}, cache, refresh)
    return {query: sorted(e['symbol'])[0] for (_, query), e in hits.items() if e['symbol']}


def resolve_endpoint(row, side, annotations, external):
    """-> (node ids, kind, mapping status)."""
    uri = endpoint_uri(row, side)
    db, identifier = identity(uri, row.get(side + '_id_database', ''), row.get(side + '_identifier', ''))
    if db in EXCLUDED or uri.startswith(WP_RDF_PREFIXES):
        return [], 'other', 'excluded_database'
    is_chemical, is_gene, genes = classify_wp_node(
        db, identifier, annotations.get(uri, {'genes': set(), 'types': set()}))
    if is_chemical and is_gene:
        kind, status = 'other', 'conflicting_type_annotations'
    elif is_chemical:
        kind, status = 'metabolite', 'metabolite_identifier_retained'
    elif is_gene:
        tag = 'mapped'
        if not genes and external.get((db, identifier), {}).get('ensembl'):
            genes, tag = set(external[(db, identifier)]['ensembl']), 'mapped_via_external_id_mapping'
        if genes:
            return sorted(genes), 'protein', tag if len(genes) == 1 else tag + '_multiple'
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


# ---------------------------------------------------------------- reporting

def print_checkpoints(table, audit, retained):
    """Unique endpoint -> Ensembl -> brain counts. Reporting only; nothing is filtered here."""
    direction = table.directionality
    print('\n===== INTERACTIONS -> ENSEMBL -> BRAIN EXPRESSION =====')
    print(f'[1] Retrieved records after the raw-export filters: {len(table):,}')
    print(f'    Unique interaction URIs: {table.interaction_uri.nunique():,}')
    print(f'    Directed records: {int(direction.eq("directed").sum()):,}; '
          f'incomplete: {int(direction.eq("incomplete").sum()):,}; '
          f'participant-only/unspecified: {int(direction.eq("undirected_or_unspecified").sum()):,}')

    for scope, rows in [('ALL source/target endpoints', audit),
                        ('DIRECTED records only', [a for a in audit if a['directionality'] == 'directed'])]:
        products, metabolites, others, mapped, brain_products = set(), set(), set(), set(), set()
        ensembl, by_status = set(), defaultdict(set)
        for a in rows:
            original, kind = a['original_uri'], a['node_type']
            if not original:
                continue
            if kind == 'metabolite':
                metabolites.add(a['node_id'] or original)
            elif kind != 'protein':
                others.add(original)
            else:
                products.add(original)
                gid = gene_id(a['node_id'])
                if gid:
                    mapped.add(original)
                    ensembl.add(gid)
                    by_status[a['brain_filter_status']].add(gid)
                    if a['brain_filter_status'] == 'brain_expressed':
                        brain_products.add(original)
        brain = by_status['brain_expressed']
        below, no_data = by_status['below_threshold'] - brain, by_status['no_expression_data'] - brain
        print(f'\n    Scope: {scope}')
        print(f'[2] BEFORE Ensembl mapping: {len(products):,} unique gene/protein endpoints; '
              f'{len(metabolites):,} unique metabolite identifiers; {len(others):,} other/unclassified endpoints.')
        print(f'[3] ENSEMBL MAPPING: {len(mapped):,} of {len(products):,} gene/protein endpoints mapped; '
              f'{len(products - mapped):,} unmapped.')
        print(f'    Mapping yields {len(ensembl):,} UNIQUE Ensembl gene IDs.')
        print(f'[4] BRAIN FILTER: {len(brain):,} of {len(ensembl):,} unique Ensembl gene IDs pass the threshold.')
        print(f'    Below threshold: {len(below - no_data):,}; no expression data: {len(no_data - below):,}; '
              f'other/unreported status: {len(ensembl - brain - below - no_data):,}.')
        print(f'    These genes represent {len(brain_products):,} original endpoints with at least one passing mapping.')

    ends = [(r[s + '_type'], r[s + '_id']) for r in retained for s in ('source', 'target')]
    genes = {gene_id(i) for t, i in ends if t == 'protein'} - {''}
    print(f'\n[5] FINAL retained network: {len(genes):,} unique Ensembl genes and '
          f'{len({i for t, i in ends if t == "metabolite"}):,} unique metabolite identifiers.')
    print('    Final counts can be lower: a passing gene also needs a retained interaction partner.')
    print(f'    Exclusion-list matches: {len(DELETED_NODES):,} unique endpoints.')
    print('Counting: source/target endpoints only, deduplicated; participant-only nodes are not counted. '
          'Before mapping, endpoints are counted by original URI, not by gene: several identifiers can map '
          'to one gene, and one endpoint can map to several genes.')
    print('=====================================================\n')


# ---------------------------------------------------------------- pipeline

def run(args):
    out = args.output_folder
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'export.status.json', {'complete': False, 'state': 'running'})

    brain_path = args.brain_expression or args.cache / 'rna_brain_region_hpa.tsv.zip'
    brain_genes, brain_summary, brain_details = load_brain_expression(
        brain_path, args.brain_threshold, args.refresh_brain_expression, args.brain_expression is None)
    measured = set(brain_summary.loc[brain_summary.measured_regions.gt(0), 'ensembl_id'])
    brain_summary.to_csv(out / 'brain_expression_gene_summary.csv', index=False)
    write_json(out / 'brain_expression_provenance.json', brain_details, indent=2)

    raw_path = out / 'wikipathways_raw_interactions.csv'
    if export_raw(args.endpoint, raw_path, args.cache / 'pathways', args.page_size, args.refresh,
                  args.pathway, args.exclude_metabolites):
        raise RuntimeError('WikiPathways retrieval was incomplete. See the raw export status and rerun without --refresh.')
    table = pd.read_csv(raw_path, dtype=str, keep_default_na=False)
    if not table.organism.str.strip().eq('Homo sapiens').all():
        raise ValueError('Non-human or unspecified organism found in raw export')
    records = table.to_dict('records')

    uris = {endpoint_uri(r, s) for r in records for s in ('source', 'target')} - {''}
    annotations = fetch_annotations(uris, args.cache / 'annotations', args.endpoint, args.refresh)
    external = map_wp_genes_to_ensembl(pending_gene_lookups(table, annotations),
                                       args.cache / 'mygene', args.refresh)

    retained, rejected, audit = [], [], []
    ok = {'brain_expressed', 'metabolite_exempt'}
    for number, row in enumerate(records, 2):
        ends, kinds = [], []
        for side in ('source', 'target'):
            ids, kind, mapping = resolve_endpoint(row, side, annotations, external)
            ends.append(ids)
            kinds.append(kind)
            for node in ids or ['']:
                audit.append({'input_csv_row': number, 'side': side, 'directionality': row['directionality'],
                              'original_uri': endpoint_uri(row, side), 'original_label': row.get(side + '_label', ''),
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
                             'source_id': a, 'target_id': b, 'source_type': kinds[0], 'target_type': kinds[1],
                             'source_ensembl': gene_id(a), 'target_ensembl': gene_id(b),
                             'source_brain_status': sa, 'target_brain_status': sb,
                             'interaction_category': interaction_category(*kinds)})

    save_csv(out / 'endpoint_audit.csv', audit,
             ['input_csv_row', 'side', 'original_uri', 'original_label', 'node_id',
              'node_type', 'mapping_status', 'brain_filter_status'])
    save_csv(out / 'removed_interactions.csv', rejected,
             ['input_csv_row', 'source_id', 'target_id', 'reason', 'interaction_uri'])
    print_checkpoints(table, audit, retained)

    genes = {r[s + '_ensembl'] for r in retained for s in ('source', 'target')} - {''}
    symbols = ensembl_symbols(genes, args.cache / 'mygene', args.refresh)
    for r in retained:
        for side in ('source', 'target'):
            if r[side + '_type'] == 'protein':  # original labels may span several genes
                r[side + '_label'] = symbols.get(r[side + '_id'], r[side + '_id'])
    sides = ('source', 'target')
    save_csv(out / 'wp_brain_interactions.csv', retained, COLUMNS + [
        f'{s}_{f}' for s in sides for f in ('id', 'type', 'ensembl', 'original_label')] +
        ['source_brain_status', 'target_brain_status', 'interaction_category'])

    # Undirected unique-pair table for network tools.
    nodes, pairs = {}, {}
    for r in retained:
        for side in sides:
            node = nodes.setdefault(r[side + '_id'], {'labels': set(), 'type': r[side + '_type']})
            if node['type'] != r[side + '_type']:
                raise ValueError(f'Conflicting endpoint types for {r[side + "_id"]}')
            node['labels'].add(r[side + '_label'])
        pair = pairs.setdefault(tuple(sorted((r['source_id'], r['target_id']))), {'pathways': set(), 'types': set()})
        pair['pathways'].update(split(r['pathway_id']))
        pair['types'].update(split(r['interaction_type']))
    edges = []
    for (a, b), item in sorted(pairs.items()):
        edges.append({'source_id': a, 'source_label': joined(nodes[a]['labels']),
                      'source_type': nodes[a]['type'], 'source_ensembl': gene_id(a),
                      'target_id': b, 'target_label': joined(nodes[b]['labels']),
                      'target_type': nodes[b]['type'], 'target_ensembl': gene_id(b),
                      'interaction_category': interaction_category(nodes[a]['type'], nodes[b]['type']),
                      'database_count': 1, 'source_database': 'WikiPathways',
                      'pathway_id': joined(item['pathways']), 'pathway_count': len(item['pathways']),
                      'interaction_types': joined(item['types'])})
    save_csv(out / 'wp_brain_edges.csv', edges,
             ['source_id', 'source_label', 'source_type', 'source_ensembl',
              'target_id', 'target_label', 'target_type', 'target_ensembl',
              'interaction_category', 'database_count', 'source_database',
              'pathway_id', 'pathway_count', 'interaction_types'])

    write_json(out / 'export.status.json', {
        'complete': True, 'organism': 'Homo sapiens', 'raw_records': len(records),
        'retained_interaction_records': len(retained), 'unique_undirected_edges': len(edges),
        'removed_pair_records': len(rejected), 'brain_threshold_nTPM': args.brain_threshold,
        'brain_expressed_genes': len(brain_genes), 'metabolites_exempt': True,
        'finished_utc': datetime.now(timezone.utc).isoformat(),
        'note': 'wp_brain_interactions.csv preserves direction; wp_brain_edges.csv is undirected. '
                'Multiple gene mappings are expanded and filtered individually; see endpoint_audit.csv.'}, indent=2)
    print(f'Finished: {len(retained):,} interaction records; {len(edges):,} unique undirected edges. {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-folder', type=Path, default=Path('wp_brain_output'))
    parser.add_argument('--cache', type=Path, default=Path('wp_brain_cache'))
    parser.add_argument('--endpoint', default='https://sparql.wikipathways.org/sparql')
    parser.add_argument('--pathway', action='append', help='Optional WP ID; repeatable')
    parser.add_argument('--page-size', type=int, default=2000)
    parser.add_argument('--refresh', action='store_true', help='Refresh pathway and identifier caches')
    parser.add_argument('--exclude-metabolites', type=Path,
                        default=Path(__file__).with_name('metabolites to delete.xlsx'),
                        help='Metabolite list CSV/XLSX (default: "metabolites to delete.xlsx" beside the script)')
    parser.add_argument('--brain-expression', type=Path,
                        help='Local HPA regional TSV/TSV.ZIP; otherwise the table is downloaded')
    parser.add_argument('--brain-threshold', type=float, default=1.0)
    parser.add_argument('--refresh-brain-expression', action='store_true')
    args = parser.parse_args()
    if args.page_size < 1:
        parser.error('--page-size must be positive')
    try:
        run(args)
    except Exception as exc:
        args.output_folder.mkdir(parents=True, exist_ok=True)
        write_json(args.output_folder / 'export.status.json', {
            'complete': False, 'error': f'{type(exc).__name__}: {exc}',
            'note': 'Do not use outputs from this folder until a run reports complete=true.'}, indent=2)
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
