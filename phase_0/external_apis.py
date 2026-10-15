import json
import urllib.request
import urllib.error
import ssl

ctx = ssl.create_default_context()

def fetch_json(url):
    req = urllib.request.Request(
        url, 
        headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"}
    )
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=20) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"error": str(e)}

def head_url(url):
    req = urllib.request.Request(
        url, 
        headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"},
        method="HEAD"
    )
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=20) as resp:
            return dict(resp.headers)
    except Exception as e:
        return {"error": str(e)}

results = {}

print("1. Fetching Figshare 27269946...")
results["figshare_27269946"] = fetch_json("https://api.figshare.com/v2/articles/27269946")

print("2. Fetching CELLxGENE 071b706a-7ea7-47a4-bddf-6457725839fc...")
results["cellxgene_071b706a"] = fetch_json("https://api.cellxgene.cziscience.com/curation/v1/collections/071b706a-7ea7-47a4-bddf-6457725839fc")

print("3. Fetching CELLxGENE f6c50495-3361-40ed-a819-fb9644396ed9...")
results["cellxgene_f6c50495"] = fetch_json("https://api.cellxgene.cziscience.com/curation/v1/collections/f6c50495-3361-40ed-a819-fb9644396ed9")

print("4. Fetching NCBI GEO GSE318789...")
geo_search = fetch_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE318789[ACCN]&retmode=json")
results["geo_search"] = geo_search
if "esearchresult" in geo_search and geo_search["esearchresult"].get("idlist"):
    geo_id = geo_search["esearchresult"]["idlist"][0]
    results["geo_summary"] = fetch_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={geo_id}&retmode=json")

print("5. Fetching PMC13318546...")
pmc_summary = fetch_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pmc&id=13318546&retmode=json")
results["pmc_summary"] = pmc_summary

with open("external_recon_results.json", "w") as f:
    json.dump(results, f, indent=2)

print("Saved external_recon_results.json successfully!")
