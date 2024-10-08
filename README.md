# ProtDomRetriever

ProtDomRetriever is a simple Python tool for retrieving protein domain information from the InterPro database based on UniProtKB accessions and specified InterPro entries. It parses InterPro JSON data and retrieves domain positions from a protein dataset.

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
   git clone https://github.com/yourusername/ProtDomRetriever.git
2. Navigate to the project directory:
   cd ProtDomRetriever
3. Install required packages:
   pip install -r requirements.txt

## Usage

Run the script using Python:

python ProtDomRetriever.py

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

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have any questions, please open an issue on the GitHub repository.

## Author

Nicolas-Frédéric Lipp, PhD  
https://github.com/NicoFrL

## Acknowledgements

This project was developed with the assistance of AI language models, which provided guidance on code structure, best practices, and documentation. The core algorithm and scientific approach were designed and implemented by the author.

