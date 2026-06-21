import os
import urllib.request

def download_file(url, filepath):
    print(f"Downloading {url} -> {filepath}")
    try:
        # Include User-Agent header to avoid block
        req = urllib.request.Request(
            url, 
            headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
        )
        with urllib.request.urlopen(req) as response:
            with open(filepath, 'wb') as f:
                f.write(response.read())
        print("Success.")
    except Exception as e:
        print(f"Failed to download: {e}")

def main():
    data_dir = r"d:\Parallel and distribution program\data"
    os.makedirs(data_dir, exist_ok=True)
    
    # NCBI Entrez efetch URLs
    genomes = {
        "sars_cov_2.fasta": "NC_045512.2",
        "sars_cov_1.fasta": "NC_004718.3",
        "mers_cov.fasta": "NC_019843.3"
    }
    
    for filename, accession in genomes.items():
        filepath = os.path.join(data_dir, filename)
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
        download_file(url, filepath)

if __name__ == "__main__":
    main()
