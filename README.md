# ProtDomRetriever

ProtDomRetriever is a simple Python tool for retrieving protein domain information from the InterPro database based on UniProtKB accessions and specified InterPro entries. The script utilizes the application programming interface (API) of InterPro, extracts the position of every domain for each entry and selects the longest domain if multiple entries overlap. Facultatively the program returns a trimmed fasta file imported from UniProt. The program allows the retrieval of multiple domains in tandem if any, and it attributes a domain number to the uniprot accession code.

Created by Nicolas-Frédéric Lipp, PhD.


## Features

- Retrieve domain information for multiple UniProtKB accessions
- Filter domains based on specified InterPro entries
- Generate TSV output with domain ranges
- Create FASTA files for the retrieved protein domains
- User-friendly GUI for file selection

## Requirements

- Python 3.6+
- Required Python packages:
  - tkinter
  - requests

## Installation

1. Clone this repository:
   git clone https://github.com/NicoFrL/ProtDomRetriever.git
2. Navigate to the project directory:
   cd ProtDomRetriever
3. Install required packages:
   pip install -r requirements.txt

## Usage

Run the script using Python:

python3 ProtDomRetriever.py

Follow the on-screen prompts to:
1. Select an input file containing UniProtKB accessions
2. Enter InterPro entries for domain filtering
3. Choose whether to fetch FASTA files for the protein domains

## Output

The script generates three main output files in a new directory:

1. `*_result_table.tsv`: A tab-separated file containing protein accessions, InterPro entries, and domain ranges
2. `*_domain_ranges.txt`: A text file listing the domain ranges for each protein
3. `*_output_domains.fasta`: A FASTA file containing the sequences of the retrieved protein domains (if FASTA retrieval is selected)

## Examples

Two example datasets are provided in the `examples` directory:

1. ORP dataset (`example1`)
2. Spectrin dataset (`example2`)

Each example includes input files, suggested InterPro entries, and sample output files.

## Manual Sequence Retrieval

Users can always use the content of the output file `*_domain_ranges.txt` at https://www.uniprot.org/id-mapping to map UniProtKB AC/ID to UniProtKB and retrieve the sequences manually, for instance as a comprehensive Excel file.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

## License

This project is distributed under **a Custom Academic and Non-Commercial License**.  
It is free to use for educational, research, and non-profit purposes.  
For commercial use, please refer to the [LICENSE](./LICENSE) file or contact the author for more information.

## Support

If you encounter any problems or have any questions, please open an issue on the GitHub repository.

## Author

Nicolas-Frédéric Lipp, PhD  
https://github.com/NicoFrL

## Development Notes
This project was developed with the assistance of AI language models, which provided guidance on code structure, best practices, and documentation. The core algorithm and scientific approach were designed and implemented by the author on the basis of InterPro documentation.

