import tkinter as tk
from tkinter import filedialog
import json
import os
import sys
import csv
import requests
import time
import logging
from queue import Queue
from threading import Thread

# Logging Setup:
# - Configure logging to display messages with timestamp and severity level
# - Show all logs of level INFO and above (INFO, WARNING, ERROR, CRITICAL)
# - Will output to console/stderr by default
# - Create a logger instance named after this module for better traceability
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Maximum number of parallel requests allowed by the InterPro API
# This limit is set to comply with InterPro's usage guidelines and prevent overloading their servers
MAX_PARALLEL_REQUESTS = 25
# Initial delay between API requests in seconds
# This helps to avoid hitting InterPro's rate limits by spacing out our requests
INITIAL_DELAY = 0.1
# Maximum number of retry attempts for failed requests
# This allows the program to handle temporary network issues or brief API unavailability
MAX_RETRIES = 3

def select_file(title, filetypes):
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title=title,
        filetypes=filetypes
    )
    return os.path.abspath(file_path) if file_path else None

def get_interpro_entries():
    entries = []
    print("\nEnter InterPro entries. You can enter entries from:")
    print("InterPro (IPR), CATH-Gene3D, CDD, HAMAP, PANTHER, Pfam, PIRSF, PRINTS,")
    print("PROSITE Patterns, PROSITE Profiles, SMART, SFLD, SUPERFAMILY, NCBIfam")
    print("Enter each entry on a new line. When finished, enter 'done' or press Enter on an empty line.")
    
    while True:
        entry = input("Enter an InterPro entry (or 'done' to finish): ").strip()
        if entry.lower() == 'done' or entry == "":
            break
        entries.append(entry)
    
    return entries

def create_output_directory(input_file):
    directory = os.path.dirname(input_file)
    input_filename = os.path.splitext(os.path.basename(input_file))[0]
    output_dir_base = f"{input_filename}_output"
    output_dir = output_dir_base
    counter = 1
    while os.path.exists(os.path.join(directory, output_dir)):
        output_dir = f"{output_dir_base}{counter}"
        counter += 1
    os.makedirs(os.path.join(directory, output_dir))
    return os.path.join(directory, output_dir)

def get_output_filenames(input_file, output_dir):
    input_filename = os.path.splitext(os.path.basename(input_file))[0]
    tsv_output = os.path.join(output_dir, f"{input_filename}_result_table.tsv")
    domain_ranges_output = os.path.join(output_dir, f"{input_filename}_domain_ranges.txt")
    fasta_output = os.path.join(output_dir, f"{input_filename}_output_domains.fasta")
    return tsv_output, domain_ranges_output, fasta_output

def worker(q, interpro_entries, results, max_domains):
    while True:
        accession = q.get()
        if accession is None:
            break
        process_accession(accession, interpro_entries, results, max_domains)
        q.task_done()

def process_accession(accession, interpro_entries, results, max_domains):
    # Construct the URL for the InterPro API request
    url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{accession}?page_size=200"
    
    # Fetch data from the API
    payload = fetch_data(url)
    
    if payload:
        # Extract domain information for each InterPro entry
        domains_by_entry = extract_domains(payload, interpro_entries)
        
        # Choose the best non-overlapping domains across all entries
        # Returns a formatted string showing which entries contributed which domains
        # and the list of selected domains sorted by position
        entry_string, chosen_domains = choose_best_domains(domains_by_entry, accession)
        
        # Update the maximum number of domains found across all proteins
        # This is used to determine how many columns needed in the output
        update_max_domains(max_domains, chosen_domains)
        
        # Prepare the result row with accession, entry string (now shows domain distribution),
        # and the chosen domains
        result_row = prepare_result_row(accession, entry_string, chosen_domains)
        
        # Store the result for this accession in the results dictionary
        results[accession] = result_row
    else:
        # If data fetch failed, store a placeholder result
        results[accession] = [accession, "N/A"]

def fetch_data(url):
    delay = INITIAL_DELAY
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, timeout=300)  # 5-minute timeout
            response.raise_for_status()
            if response.status_code == 408:  # Request Timeout
                print(f"Request timeout for {url}. Consider breaking into smaller queries.")
                return None
            return response.json()
        except requests.RequestException as e:
            print(f"Error fetching data (attempt {attempt + 1}): {e}")
            if attempt < MAX_RETRIES - 1:
                delay *= 2  # Exponential backoff
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print(f"Failed to fetch data after {MAX_RETRIES} attempts")
                return None
        finally:
            time.sleep(INITIAL_DELAY)  # Always wait a bit between requests

def extract_domains(payload, interpro_entries):
    # Create a dictionary to store domains for each InterPro entry
    domains_by_entry = {entry: [] for entry in interpro_entries}
    
    # Iterate through the API results
    for item in payload.get("results", []):
        if item["metadata"]["accession"] in interpro_entries:
            protein_data = item["proteins"][0]
            locations = protein_data.get("entry_protein_locations", [])
            
            # Extract domain ranges from the location data
            for location in locations:
                fragments = location.get("fragments", [])
                for fragment in fragments:
                    start = fragment.get("start", "N/A")
                    end = fragment.get("end", "N/A")
                    if start != "N/A" and end != "N/A":
                        domains_by_entry[item["metadata"]["accession"]].append((int(start), int(end)))
    
    return domains_by_entry

def choose_best_domains(domains_by_entry, accession):  # Add accession parameter:
    """
    Select the best domains across all entries, taking the longest when there's overlap
    """
    # Collect all domains from all entries with their source
    all_domains = []
    entries_by_domain = {}  # To track which entry each domain came from
    
    for entry, domains in domains_by_entry.items():
        for domain in domains:
            all_domains.append(domain)
            entries_by_domain[domain] = entry

    # Sort domains by length (longest first)
    all_domains.sort(key=lambda x: x[1] - x[0], reverse=True)
    
    # Final selected domains and their entries
    selected_domains = []
    entry_domain_map = {}  # Map entries to their domain numbers
    domain_counter = 1  # To give each domain a unique number
    
    def domains_overlap(d1, d2):
        return not (d1[1] < d2[0] or d1[0] > d2[1])

    for domain in all_domains:
        # Check if this domain overlaps with any selected domain
        overlap = False
        for selected in selected_domains:
            if domains_overlap(domain, selected):
                overlap = True
                break
        
        if not overlap:
            selected_domains.append(domain)
            entry = entries_by_domain[domain]
            if entry not in entry_domain_map:
                entry_domain_map[entry] = []
            entry_domain_map[entry].append(f"d{domain_counter}")
            domain_counter += 1

    # Format the entry string with domain numbers
    entry_parts = []
    for entry, domains in entry_domain_map.items():
        ranges = []
        for domain in domains:
            # Get the index of this domain (remove 'd' and convert to int - 1)
            idx = int(domain[1:]) - 1
            # Get the actual range
            start, end = selected_domains[idx]
            ranges.append(f"{domain}:[{start}, {end}]")
        entry_parts.append(f"{entry} ({','.join(ranges)})")
    
    entry_string = " + ".join(entry_parts)
    
    # Include the accession number in the log message
    if selected_domains:  # Only add accession if we found domains
        logger.info(f"Selected domains from entries: {entry_string} in {accession}")
    else:
        logger.info(f"No domains found for {accession}")
    
    return entry_string, selected_domains

def update_max_domains(max_domains, chosen_domains):
    # Update the maximum number of domains found across all proteins
    # max_domains[0] is used because it's passed as a mutable list to allow updating across threads
    max_domains[0] = max(max_domains[0], len(chosen_domains))

def prepare_result_row(accession, chosen_entry, chosen_domains):
    # Initialize the result row with the protein accession and the chosen InterPro entry
    result_row = [accession, chosen_entry]
    
    # Add the start and end positions of each domain to the result row
    for domain in chosen_domains:
        result_row.extend(domain) # domain is a tuple of (start, end), so we extend the list with both values
    # Return the complete result row
    return result_row

def submit_id_mapping(ids):
    url = "https://rest.uniprot.org/idmapping/run"
    data = {
        'ids': ids,
        'from': "UniProtKB_AC-ID",
        'to': "UniProtKB"
    }
    response = requests.post(url, data=data)
    return response.json()['jobId']

def check_job_status(job_id, timeout=300):
    url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    start_time = time.time()
    while True:
        response = requests.get(url)
        status = response.json()
        
        if 'jobStatus' in status:
            if status['jobStatus'] == "FINISHED":
                return status
        elif 'results' in status or 'failedIds' in status:
            return status

        if time.time() - start_time > timeout:
            print(f"Job timed out after {timeout} seconds. Please check the status manually:")
            print(f"https://www.uniprot.org/id-mapping/uniprotkb/{job_id}/overview")
            return None
        
        current_time = time.time()
        if current_time - start_time >= 15:
            print("Still retrieving data... Please wait.")
            start_time = current_time
        time.sleep(1)

def get_job_results(job_id):
    details_url = f"https://rest.uniprot.org/idmapping/details/{job_id}"
    response = requests.get(details_url)
    details = response.json()
    # Error handling: Log and display error if redirect URL is missing from API response
    if 'redirectURL' not in details:
        logger.error("Error: Couldn't get job results. Please check manually:")
        logger.error(f"https://www.uniprot.org/id-mapping/uniprotkb/{job_id}/overview")
        return None
    
    results_url = details['redirectURL']
    params = {
        'format': 'fasta',
        'size': '500'
    }
    # Initialize variable to store all results from potentially multiple pages
    all_results = ""
    while True:
        # Make HTTP request and handle potential errors
        response = requests.get(results_url, params=params)
        if response.status_code != 200:
            logger.error(f"Error fetching results: HTTP {response.status_code}")
            return None
        # Get content from current page
        page_content = response.text
        if not page_content.strip():
            break  # No more results
        
        all_results += page_content
        
        # Check if there are more results (Check for pagination: are there more pages?)
        if 'Link' in response.headers:
            next_link = [link.strip() for link in response.headers['Link'].split(',') if 'next' in link]
            if next_link:
                results_url = next_link[0].split(';')[0].strip('<>')
            else:
                break
        else:
            break
    # Log appropriate messages about the results
    if not all_results:
        logger.warning("No results were fetched. The output file may be empty.")
    else:
        logger.info("Successfully fetched all results")
    
    return all_results

def parse_fasta(fasta_content, domain_ranges):
    sequences = {}
    current_accession = None
    current_sequence = ""

    for line in fasta_content.split('\n'):
        if line.startswith('>'):
            if current_accession:
                sequences[current_accession] = current_sequence
            current_accession = line.split('|')[1]
            current_sequence = ""
        else:
            current_sequence += line

    if current_accession:
        sequences[current_accession] = current_sequence

    logger.info(f"Parsed {len(sequences)} sequences from FASTA content") # Log how many sequences are found

    domain_sequences = {}
    for accession, ranges in domain_ranges.items():
        if accession in sequences:
            for start, end in ranges:
                try:
                    domain_key = f"{accession}[{start}-{end}]" # Create the domain key (like "P12345[10-50]")
                    domain_sequences[domain_key] = sequences[accession][int(start)-1:int(end)] # Try to extract the sequence
                    logger.info(f"Extracted domain {domain_key}") # Log success
                except IndexError:
                    logger.warning(f"Failed to extract domain {accession}[{start}-{end}]: Index out of range") # Log when the start/end positions are invalid
        else:
            logger.warning(f"Accession {accession} not found in FASTA content")  # Log when the accession number isn't found

    logger.info(f"Extracted {len(domain_sequences)} domain sequences")
    return domain_sequences

def main():
    # Select input file
    input_file = select_file("Select input Accession List file", [("Text files", "*.txt")])
    if not input_file:
        print("No file selected. Exiting.")
        return

    print(f"Selected file: {input_file}")

    # Get InterPro entries from user
    interpro_entries = get_interpro_entries()

    if not interpro_entries:
        print("No InterPro entries provided. Exiting.")
        return

    print(f"\nProcessing with the following InterPro entries: {', '.join(interpro_entries)}")

    # Create output directory and get output filenames
    output_dir = create_output_directory(input_file)
    tsv_output, domain_ranges_output, fasta_output = get_output_filenames(input_file, output_dir)
    print(f"Results will be saved to: {output_dir}")

    # Read accessions from input file
    with open(input_file, "r") as file:
        accessions = [line.strip() for line in file if line.strip()]

    results = {}
    max_domains = [0]  # Use a list instead of Value

    # Set up multi-threading for parallel processing
    q = Queue()
    threads = []
    for _ in range(MAX_PARALLEL_REQUESTS):
        t = Thread(target=worker, args=(q, interpro_entries, results, max_domains))
        t.start()
        threads.append(t)

    # Add accessions to the queue for processing
    for accession in accessions:
        q.put(accession)

    # Wait for all tasks to complete
    q.join()

    # Stop worker threads
    for _ in range(MAX_PARALLEL_REQUESTS):
        q.put(None)
    for t in threads:
        t.join()

    # Prepare headers for output files
    header_parts = ["Protein Accession", "InterPro Entry"]
    for i in range(1, max(max_domains[0], 1) + 1):
        header_parts.extend([f"Start {i}", f"End {i}"])
    
    header_tsv = header_parts

    # Write results to output files
    with open(tsv_output, "w", newline='') as tsv_file, \
         open(domain_ranges_output, "w") as domain_ranges_file:
        
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        
        tsv_writer.writerow(header_tsv)

        domain_ranges = {}
        for accession in accessions:
            result = results.get(accession, [accession, "N/A", "N/A", "N/A"])
            if len(result) > 4:  # If domains were found
                row = result + ["N/A"] * (2 * max_domains[0] - len(result) + 2)
            else:  # If no domains were found
                row = result
            
            tsv_writer.writerow(row)

            # Write to domain ranges file and update domain_ranges
            if len(result) > 2:
                domain_ranges[accession] = []
                for i in range(2, len(result), 2):
                    start, end = result[i], result[i+1]
                    if start != "N/A" and end != "N/A":
                        domain_ranges_file.write(f"{accession}[{start}-{end}]\n")
                        domain_ranges[accession].append((int(start), int(end)))

    print(f"\nProcessing complete. Results saved to {output_dir}")

    # FASTA retrieval
    fetch_fasta = input("Do you want to fetch FASTA files for the protein domains? (Y/N): ").upper()
    if fetch_fasta == 'Y':
        job_id = submit_id_mapping(','.join(domain_ranges.keys()))
        print(f"Job submitted with ID: {job_id}")
        print(f"URL: https://www.uniprot.org/id-mapping/uniprotkb/{job_id}/overview")

        print("Retrieving data... This may take a while.")
        status = check_job_status(job_id)
        if status is None:
            return

        print("Downloading FASTA data...")
        fasta_content = get_job_results(job_id)
        if fasta_content is None:
            return

        domain_sequences = parse_fasta(fasta_content, domain_ranges)

        with open(fasta_output, "w") as fasta_file:
            for domain, sequence in domain_sequences.items():
                fasta_file.write(f">{domain}\n{sequence}\n")

        print(f"FASTA results have been saved to {fasta_output}")

if __name__ == "__main__":
    main()
