#!/usr/bin/env python3
"""
Extract AWS S3, EBI FTP, and HTTP links for SRA/ENA accessions using Selenium.

Usage:
    ./obtainrawfastqdownloadlinks.py --accession SRR32578126
    ./obtainrawfastqdownloadlinks.py --list accessions.txt --threads 16 [--verbose]

Requirements:
    - Python 3
    - Selenium, Firefox, geckodriver, and required Conda packages.
"""

import os
import sys
import time
import argparse
import logging
import concurrent.futures
import multiprocessing
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from tqdm import tqdm

# Constants
RETRIES = 3
OUTPUT_FILE = "rawfastq_download_links.txt"

# Suppress DEBUG logs from Selenium and urllib3
logging.getLogger("selenium").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

def setup_logging(log_file="scraper.log", verbose=False):
    """
    Sets up two handlers:
      - A file handler that logs INFO and above (no DEBUG).
      - A console handler that logs INFO messages only if verbose is enabled,
        and is WARNING if not verbose.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)  # Only allow INFO+ across the board
    logger.handlers.clear()

    # File handler: logs INFO and above to scraper.log
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Console handler: show INFO if verbose, otherwise WARNING
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO if verbose else logging.WARNING)
    console_formatter = logging.Formatter("%(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

def initialize_output_file():
    """
    Create/overwrite the output file with a header.
    """
    with open(OUTPUT_FILE, "w") as f:
        f.write("accession\tdownload_path\n")

# Set up Firefox for headless execution
firefox_options = Options()
firefox_options.add_argument("--headless")
firefox_options.add_argument("--no-sandbox")
firefox_options.add_argument("--disable-dev-shm-usage")

def get_sra_links(accession):
    """
    Scrape the NCBI Run Browser page for download links for the given accession.
    Returns a list of valid links, or None if none are found.
    """
    ncbi_url = f"https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc={accession}&display=data-access"
    attempt = 1
    sra_links = None

    while attempt <= RETRIES:
        driver = None
        try:
            # Removed debug logs from here – they won't appear in the file anyway
            service = Service()
            driver = webdriver.Firefox(service=service, options=firefox_options)
            driver.get(ncbi_url)

            # Wait for table to load
            WebDriverWait(driver, 60).until(
                EC.presence_of_element_located(
                    (By.XPATH, '//*[@id="ph-run-browser-data-access"]/div[2]/table')
                )
            )
            time.sleep(5)  # Extra time for JavaScript to render

            # Find rows in the table
            rows = driver.find_elements(By.XPATH, '/html/body/div[1]/div[4]/div/div[3]/div[4]/div[2]/table/tbody/tr')
            if not rows:
                sra_links = None
            else:
                sra_links = []
                for row in rows:
                    try:
                        link_element = row.find_element(By.XPATH, "./td[6]")
                        link_text = link_element.text.strip()
                        # If there's an anchor, extract href
                        try:
                            a_tag = row.find_element(By.XPATH, "./td[6]/a")
                            link_text = a_tag.get_attribute("href")
                        except Exception:
                            pass
                        if link_text.startswith(("s3://", "ftp://", "http://", "https://")):
                            sra_links.append(link_text)
                    except Exception:
                        pass

            if sra_links:
                return sra_links
        except Exception:
            pass
        finally:
            if driver:
                driver.quit()

        attempt += 1
        if attempt <= RETRIES:
            time.sleep(5)
    return sra_links

def save_links(accession, links):
    """
    Append the accession and each link to the output file.
    """
    try:
        with open(OUTPUT_FILE, "a") as f:
            for link in links:
                f.write(f"{accession}\t{link}\n")
    except Exception as e:
        logging.error(f"Error saving links for {accession}: {e}")

def process_accession(accession):
    """
    Process a single accession. Returns (accession, True/False).
    """
    links = get_sra_links(accession)
    if links:
        save_links(accession, links)
        return (accession, True)
    else:
        return (accession, False)

def get_optimal_worker_count(user_threads=None):
    """
    Determine number of worker processes: either user-specified or 75% of available cores.
    """
    available_cores = multiprocessing.cpu_count()
    default_workers = max(1, int(available_cores * 0.75))
    return user_threads if user_threads else default_workers

def main():
    parser = argparse.ArgumentParser(description="Parallelized SRA/ENA scraper for download links.")
    parser.add_argument("--accession", type=str, help="Single SRA/ENA accession.")
    parser.add_argument("--list", type=str, help="File with a list of accessions.")
    parser.add_argument("--threads", type=int, help="Number of parallel workers (default: 75%% of cores).")
    parser.add_argument("--verbose", action="store_true", help="Enable minimal console output (SUCCESS/FAILED messages).")
    args = parser.parse_args()

    setup_logging(verbose=args.verbose)
    initialize_output_file()

    accessions = []
    if args.accession:
        accessions.append(args.accession.strip())
    if args.list:
        if not os.path.exists(args.list):
            logging.error(f"List file {args.list} not found.")
            sys.exit(1)
        with open(args.list, "r") as f:
            accessions = [line.strip() for line in f if line.strip()]

    if not accessions:
        logging.error("No accession provided. Use --accession or --list.")
        sys.exit(1)

    worker_count = get_optimal_worker_count(args.threads)
    logging.info(f"Processing {len(accessions)} accessions with {worker_count} workers.")

    # In non-verbose mode, use a tqdm progress bar
    if not args.verbose:
        with tqdm(total=len(accessions), desc="Processing Accessions", dynamic_ncols=True) as progress:
            with concurrent.futures.ProcessPoolExecutor(max_workers=worker_count) as executor:
                futures = {executor.submit(process_accession, acc): acc for acc in accessions}
                for future in concurrent.futures.as_completed(futures):
                    # We don't need the result here for non-verbose
                    progress.update(1)
    else:
        # In verbose mode, print minimal messages for each accession
        with concurrent.futures.ProcessPoolExecutor(max_workers=worker_count) as executor:
            for acc, success in executor.map(process_accession, accessions):
                if success:
                    print(f"SUCCESS: {acc}")
                else:
                    print(f"FAILED: {acc}")

    logging.info("All accessions processed.")

if __name__ == "__main__":
    main()
