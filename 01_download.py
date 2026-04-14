import urllib.request
import os

os.makedirs("pdb_files", exist_ok=True)

pdbs = ['3dfr','4dfr','2hhb','1mbd','2lzm',
        '1lyz','3ptn','2gch','3cyt','3c2c','2aza','1pcy']

for pdb_id in pdbs:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out = f"pdb_files/{pdb_id}.pdb"
    print(f"Downloading {pdb_id}...", end=' ')
    urllib.request.urlretrieve(url, out)
    print("done")

print("\nAll files downloaded!")