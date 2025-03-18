#!/usr/bin/env python3
"""
obtainrawfastqdownloadlinks.py - Extract AWS S3, EBI FTP, and HTTP links for SRA/ENA accessions.

Usage:
    ./obtainrawfastqdownloadlinks.py --accession SRR28962973
    ./obtainrawfastqdownloadlinks.py --list accessions.txt

Options:
    --accession <str>   Extract links for a single SRA/ENA accession.
    --list <file>       Provide a file with a list of accessions.

Requirements:
    - Python 3
    - Conda environment with `selenium requests beautifulsoup4 firefox geckodriver`
"""

import os
import sys
import argparse
import time
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from tqdm import tqdm

# Directories
WORKDIR = os.getcwd()
DOWNLOAD_LINKS_FILE = os.path.join(WORKDIR, "download_links.txt")
LOG_FILE = os.path.join(WORKDIR, "scraper.log")

# Set up Firefox options
firefox_options = Options()
firefox_options.add_argument("--headless")  # Run in headless mode (no UI)

# Function to log messages
def log_message(accession, message):
    """Log success or failure messages."""
    with open(LOG_FILE, "a") as log:
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {accession} - {message}\n")
    print(message)

# Function to scrape all possible download links
def get_sra_links(accession):
    """Scrape NCBI Run Browser for AWS, FTP, and HTTP links using Selenium."""
    ncbi_url = f"https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc={accession}&display=data-access"

    # Initialize WebDriver (Firefox via Conda)
    service = Service()
    driver = webdriver.Firefox(service=service, options=firefox_options)
    driver.get(ncbi_url)

    # Wait for the data access table to load
    try:
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.XPATH, '//*[@id="ph-run-browser-data-access"]/div[2]/table'))
        )
    except Exception as e:
        log_message(accession, f"⚠️ Page did not load properly: {e}")
        driver.quit()
        return None

    # Extract links (AWS S3, EBI FTP, HTTP)
    sra_links = []
    for row in driver.find_elements(By.XPATH, '//*[@id="ph-run-browser-data-access"]/div[2]/table/tbody/tr'):
        try:
            # Extract download links (last column in the table)
            link_element = row.find_elements(By.TAG_NAME, "td")[-1]
            link_text = link_element.text.strip()

            # Ensure it's a valid link
            if link_text.startswith(("s3://", "ftp://", "http://", "https://")):
                sra_links.append(link_text)
        except Exception:
            pass  # Ignore if extraction fails

    driver.quit()

    if not sra_links:
        log_message(accession, "⚠️ No download links found.")
        return None

    return sra_links

# Function to save links
def save_sra_links(accession, links):
    """Save SRA links to file."""
    with open(DOWNLOAD_LINKS_FILE, "a") as file:
        for url in links:
            file.write(f"{accession}\t{url}\n")

    log_message(accession, f"✅ Download links saved.")

# Main function
def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Fetch AWS S3, EBI FTP, and HTTP links for SRA/ENA accessions.")
    parser.add_argument("--accession", type=str, help="SRA/ENA run accession.")
    parser.add_argument("--list", type=str, help="File with a list of accessions.")
    args = parser.parse_args()

    accessions = []

    if args.accession:
        accessions.append(args.accession.strip())

    if args.list:
        if not os.path.exists(args.list):
            print(f"❌ ERROR: File {args.list} not found.")
            sys.exit(1)

        with open(args.list, "r") as file:
            accessions.extend([line.strip() for line in file])

    if not accessions:
        print("❌ ERROR: Provide either --accession or --list. Use --help for usage.")
        sys.exit(1)

    for accession in tqdm(accessions, desc="Processing accessions"):
        log_message(accession, f"🚀 Fetching download links for {accession}")

        links = get_sra_links(accession)
        if not links:
            log_message(accession, f"❌ No FASTQ files found for {accession}.")
            continue

        save_sra_links(accession, links)

    log_message("ALL", "🎉 Download link extraction completed.")

if __name__ == "__main__":
    main()
