import json
from kaggle.api.kaggle_api_extended import KaggleApi

api = KaggleApi()
api.authenticate()

queries = [
    "BoneMarrowMap",
    "Zeng AML",
    "scAML",
    "AML scAtlas",
    "748,679",
    "748679",
    "van Galen AML",
    "van Galen",
    "GSE318789",
    "AML CRISPRi lncRNA",
    "AML CRISPRi",
    "CRISPRi lncRNA",
    "MOLM-13",
    "AML single cell"
]

results_by_query = {}
all_datasets = {}

for q in queries:
    datasets = api.dataset_list(search=q)
    results_by_query[q] = []
    for d in datasets:
        ref = f"{d.ref}"
        entry = {
            "ref": ref,
            "title": d.title,
            "size": getattr(d, 'size', 'N/A'),
            "totalBytes": getattr(d, 'totalBytes', None),
            "lastUpdated": str(getattr(d, 'lastUpdated', 'N/A')),
            "downloadCount": getattr(d, 'downloadCount', 0),
            "url": f"https://www.kaggle.com/datasets/{ref}"
        }
        results_by_query[q].append(entry)
        if ref not in all_datasets:
            all_datasets[ref] = entry

print(f"Total unique datasets found across all queries: {len(all_datasets)}")
for q, res in results_by_query.items():
    print(f"Query '{q}': {len(res)} results")
    for r in res[:5]:
        print(f"  - {r['ref']} ({r['size']}): {r['title']}")

# Save detailed results
with open("kaggle_search_results.json", "w") as f:
    json.dump({"by_query": results_by_query, "all_unique": all_datasets}, f, indent=2)
