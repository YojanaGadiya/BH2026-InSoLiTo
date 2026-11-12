"""
Load the tool pairs of registered workflows into Neo4j as WORKFLOW_NEXT relationships.

Goal: show in InSoLiTo, as their own layer, the tool pairs that registered
workflows (WorkflowHub: Galaxy, Nextflow, nf-core) run one after the other,
next to the co-citation edges (METAOCCUR) and the literature order (CITED_BEFORE).

- One relationship per row of combined_edges.csv: (source)-[:WORKFLOW_NEXT]->(target),
  with the number of workflows per source and their total. The literature columns
  of that file are not loaded (that information is already in CITED_BEFORE).
- The Nextflow column excludes nf-core workflows, so the three counts can be added.
- Tools are matched to existing Tool/Database nodes by `label` (bio.tools ID),
  ignoring case. No data is changed: only the lookup ignores case.
- Rows whose source or target is not in Neo4j are skipped (and counted).
- Only MERGE: nodes are never created or deleted. Re-running updates the
  numbers of edges already loaded (it does not remove edges that are no longer
  in the CSV).
- Only against localhost.

With DRY_RUN = True it only prints what it would load and writes nothing.

Undo (Neo4j Browser):
    :auto MATCH ()-[r:WORKFLOW_NEXT]->() CALL { WITH r DELETE r } IN TRANSACTIONS OF 10000 ROWS

Run from the repo root, with Neo4j up:
    .venv/bin/python workflows/scripts/load_workflow_edges_neo4j.py
"""
import duckdb
from neo4j import GraphDatabase

EDGES_CSV = 'workflows/data/biotools/combined/combined_edges.csv'
NEO4J_URI = 'bolt://localhost:7687'
DRY_RUN = True
BATCH_SIZE = 10_000

# Number of workflows per source, stored on each relationship
SOURCES = ['n_wf_galaxy', 'n_wf_nextflow', 'n_wf_nfcore']

if 'localhost' not in NEO4J_URI and '127.0.0.1' not in NEO4J_URI:
    raise SystemExit(f'Refusing to write to a non-local Neo4j: {NEO4J_URI}')

driver = GraphDatabase.driver(NEO4J_URI, auth=None)

# 1. Lookup table: lowercase label -> Neo4j node id
with driver.session() as session:
    nodes = session.run(
        'MATCH (n) WHERE n:Tool OR n:Database RETURN id(n) AS id, n.label AS label').data()

by_lower = {}
for n in nodes:
    by_lower.setdefault(n['label'].lower(), []).append(n)

# Labels that only differ in case would point to two nodes: leave them out
clashes = {k: [n['label'] for n in v] for k, v in by_lower.items() if len(v) > 1}
print(f'Tool/Database nodes:         {len(nodes):,}')
print(f'Labels clashing by case:     {len(clashes):,} {list(clashes.values())[:5]}')
node_id = {k: v[0]['id'] for k, v in by_lower.items() if len(v) == 1}

# 2. Read the edges and map both ends to node ids
con = duckdb.connect()
rows = con.execute(
    f"SELECT source, target, {', '.join(SOURCES)} FROM read_csv_auto('{EDGES_CSV}')").fetchall()
print(f'Rows in CSV:                 {len(rows):,}')

batch_rows, seen = [], set()
missing = self_loops = duplicates = 0
for source, target, *values in rows:
    a, b = node_id.get(str(source).lower()), node_id.get(str(target).lower())
    if a is None or b is None:
        missing += 1
        continue
    if a == b:              # two CSV IDs that land on the same node
        self_loops += 1
        continue
    if (a, b) in seen:      # same node pair reached from differently-cased IDs
        duplicates += 1
        continue
    seen.add((a, b))
    props = dict(zip(SOURCES, values))
    props['n_workflows'] = sum(values)
    batch_rows.append({'a': a, 'b': b, 'props': props})

print(f'Skipped, tool not in Neo4j:  {missing:,}')
print(f'Skipped, same node:          {self_loops:,}')
print(f'Skipped, repeated pair:      {duplicates:,}')
print(f'Edges to load:               {len(batch_rows):,}')

# 3. Load in batches
if DRY_RUN:
    print('DRY_RUN = True: nothing written.')
else:
    query = (
        'UNWIND $rows AS row '
        'MATCH (a) WHERE id(a) = row.a '
        'MATCH (b) WHERE id(b) = row.b '
        'MERGE (a)-[r:WORKFLOW_NEXT]->(b) '
        'SET r += row.props')
    with driver.session() as session:
        for i in range(0, len(batch_rows), BATCH_SIZE):
            session.run(query, rows=batch_rows[i:i + BATCH_SIZE]).consume()
            print(f'  loaded {min(i + BATCH_SIZE, len(batch_rows)):,} / {len(batch_rows):,}')
        total = session.run(
            'MATCH ()-[r:WORKFLOW_NEXT]->() RETURN count(r) AS n').single()['n']
    print(f'WORKFLOW_NEXT relationships in Neo4j: {total:,}')

driver.close()
