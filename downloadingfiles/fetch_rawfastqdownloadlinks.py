#!/usr/bin/env python3
"""
Parallelized version of fetch_rawfastqdownloadlinks.py
Extract AWS S3, EBI FTP, and HTTP links for SRA/ENA accessions using Selenium.

Usage:
    ./fetch_rawfastqdownloadlinks.py --accession SRR32578126
    ./fetch_rawfastqdownloadlinks.py --list accessions.txt --threads 16 [--verbose]
"""

import os
import sys
import time
import argparse
import logging
import concurrent.futures
import multiprocessing
from tqdm import tqdm

from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Constants
RETRIES = 3
OUTPUT_FILE = "download_links.txt"
LOG_FILE = "scraper.log"

###############################################################################
# 1) A custom flush handler so each log record is written to file immediately #
###############################################################################
class FlushFileHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

# Suppress DEBUG logs from Selenium and urllib3
logging.getLogger("selenium").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

def setup_logging(verbose=False):
    """
    Configure:
    - A file handler (FlushFileHandler) that logs INFO+ to scraper.log
    - A console handler that logs INFO if verbose, else WARNING
    """
    logger = logging.getLogger()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    # Custom file handler for real-time logs
    file_handler = FlushFileHandler(LOG_FILE, mode="w")
    file_handler.setLevel(logging.INFO)
    file_fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s",
                                 datefmt="%Y-%m-%d %H:%M:%S")
    file_handler.setFormatter(file_fmt)
    logger.addHandler(file_handler)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_level = logging.INFO if verbose else logging.WARNING
    console_handler.setLevel(console_level)
    console_fmt = logging.Formatter("%(message)s")
    console_handler.setFormatter(console_fmt)
    logger.addHandler(console_handler)

def initialize_output_file():
    """
    Create or overwrite the output file with a header.
    """
    with open(OUTPUT_FILE, "w") as f:
        f.write("accession\tdownload_path\n")

# Firefox options for headless scraping
firefox_options = Options()
firefox_options.add_argument("--headless")
firefox_options.add_argument("--no-sandbox")
firefox_options.add_argument("--disable-dev-shm-usage")

def get_sra_links(accession):
    """
    Scrape the NCBI Run Browser page for download links for the given accession.
    Return a list of valid links, or None if none found.
    """
    ncbi_url = f"https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc={accession}&display=data-access"
    attempt = 1
    sra_links = None

    while attempt <= RETRIES:
        driver = None
        try:
            service = Service()
            driver = webdriver.Firefox(service=service, options=firefox_options)
            driver.get(ncbi_url)

            # Wait for table to load
            WebDriverWait(driver, 60).until(
                EC.presence_of_element_located((By.XPATH, '//*[@id="ph-run-browser-data-access"]/div[2]/table'))
            )
            time.sleep(5)  # extra time for JS rendering

            # Find rows in the table
            rows = driver.find_elements(By.XPATH, '/html/body/div[1]/div[4]/div/div[3]/div[4]/div[2]/table/tbody/tr')
            if rows:
                sra_links = []
                for row in rows:
                    try:
                        link_element = row.find_element(By.XPATH, "./td[6]")
                        link_text = link_element.text.strip()
                        # If there's a clickable link, get its href
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
    Process a single accession. Logs SUCCESS/FAILED to the log file (info-level).
    Returns (accession, bool) indicating success/failure.
    """
    links = get_sra_links(accession)
    if links:
        save_links(accession, links)
        logging.info(f"SUCCESS: {accession}")
        return (accession, True)
    else:
        logging.info(f"FAILED: {accession}")
        return (accession, False)

def get_optimal_worker_count(user_threads=None):
    """
    Determine # of worker processes: user-specified or 75% of CPU cores.
    """
    cores = multiprocessing.cpu_count()
    default_workers = max(1, int(cores * 0.75))
    return user_threads if user_threads else default_workers

#################################################################################
# Tailing logic: read the log file in real time, skipping the first line,        #
# updating the progress bar for each "SUCCESS" or "FAILED", then stopping        #
# if we see "All accessions processed."                                         #
#################################################################################
def tail_log_and_update_progress(progress):
    with open(LOG_FILE, "r") as f:
        # Skip the first line
        _ = f.readline()

        while True:
            line = f.readline()
            if not line:
                time.sleep(0.1)  # Wait for new log lines
                continue

            if "SUCCESS" in line or "FAILED" in line:
                progress.update(1)

            if "All accessions processed." in line:
                # Done reading
                break

def main():
    parser = argparse.ArgumentParser(
        description="Parallelized SRA/ENA scraper for AWS S3, EBI FTP, HTTP links."
    )
    parser.add_argument("--accession", type=str, help="Single SRA/ENA accession.")
    parser.add_argument("--list", type=str, help="File with a list of accessions.")
    parser.add_argument("--threads", type=int, help="Number of parallel workers.")
    parser.add_argument("--verbose", action="store_true", help="Print success/fail to console.")
    args = parser.parse_args()

    setup_logging(verbose=args.verbose)
    initialize_output_file()

    # Collect accessions
    accessions = []
    if args.accession:
        accessions.append(args.accession.strip())
    if args.list:
        if not os.path.exists(args.list):
            logging.error(f"List file {args.list} not found.")
            sys.exit(1)
        with open(args.list, "r") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
            accessions.extend(lines)

    if not accessions:
        logging.error("No accessions found. Use --accession or --list.")
        sys.exit(1)

    worker_count = get_optimal_worker_count(args.threads)
    logging.info(f"Processing {len(accessions)} accessions with {worker_count} workers.")

    # Launch parallel tasks
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_count) as executor:
        futures = [executor.submit(process_accession, acc) for acc in accessions]

    # If not verbose => show a progress bar that increments each time
    # a line is written to the log for a success/fail.
    if not args.verbose:
        progress = tqdm(total=len(accessions), desc="Processing Accessions", dynamic_ncols=True)
        tail_log_and_update_progress(progress)
        progress.close()
    else:
        # If verbose => each success/fail was printed to the console already
        pass

    logging.info("All accessions processed.")

if __name__ == "__main__":
    main()
