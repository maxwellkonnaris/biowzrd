#!/usr/bin/env python3
"""
obtainrawfastqdownloadlinks.py - Extract AWS S3, EBI FTP, and HTTP links for SRA/ENA accessions.

Usage:
    ./obtainrawfastqdownloadlinks.py --accession SRR28962973
    ./obtainrawfastqdownloadlinks.py --list accessions.txt

Requirements:
    - Python 3
    - Selenium, Firefox, geckodriver, and the necessary Conda environment packages.
"""

import os
import sys
import time
import argparse
import logging
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Configure logging to both file and console
def setup_logging(log_file="scraper.log"):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

# Set up Firefox options (headless mode)
firefox_options = Options()
firefox_options.add_argument("--headless")

def get_sra_links(accession):
    """
    Scrape NCBI Run Browser for AWS, FTP, and HTTP links using Selenium.
    
    Returns a list of valid download links for the given accession or None if none found.
    """
    ncbi_url = f"https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc={accession}&display=data-access"
    driver = None
    sra_links = []
    
    try:
        # Initialize WebDriver
        service = Service()
        driver = webdriver.Firefox(service=service, options=firefox_options)
        driver.get(ncbi_url)
        logging.info(f"Accessing URL: {ncbi_url}")

        # Wait for the table element to load
        try:
            WebDriverWait(driver, 20).until(
                EC.presence_of_element_located(
                    (By.XPATH, '/html/body/div[1]/div[4]/div/div[3]/div[4]/div[2]/table')
                )
            )
        except Exception as e:
            logging.warning(f"Page did not load properly for accession {accession}: {e}")
            logging.debug(driver.page_source)
            return None

        # Extra wait for JavaScript rendering (adjust if necessary)
        time.sleep(5)

        # Find table rows using full XPath
        rows = driver.find_elements(By.XPATH, '/html/body/div[1]/div[4]/div/div[3]/div[4]/div[2]/table/tbody/tr')
        if not rows:
            logging.warning(f"No rows found in table for accession {accession}.")
            logging.debug(driver.page_source)
            return None

        # Iterate over rows and extract the link from the 6th column
        for row in rows:
            try:
                # First, attempt to get direct text from the 6th column
                link_element = row.find_element(By.XPATH, "./td[6]")
                link_text = link_element.text.strip()

                # If the column contains a clickable link (<a> tag), extract the href attribute
                try:
                    a_tag = row.find_element(By.XPATH, "./td[6]/a")
                    link_text = a_tag.get_attribute("href")
                except Exception:
                    pass  # If no <a> tag, continue with text content

                # Validate the link format
                if link_text.startswith(("s3://", "ftp://", "http://", "https://")):
                    sra_links.append(link_text)
            except Exception as e:
                logging.error(f"Error extracting link from a row for accession {accession}: {e}")
                continue

        if not sra_links:
            logging.warning(f"No valid download links found for accession {accession}.")
            return None

        return sra_links

    except Exception as e:
        logging.error(f"An error occurred for accession {accession}: {e}")
        return None
    finally:
        if driver is not None:
            driver.quit()

def save_links(accession, links, output_file="download_links.txt"):
    """
    Save the extracted download links to a file.
    
    Each line is written in the format: accession<TAB>link
    """
    try:
        with open(output_file, "a") as f:
            for link in links:
                f.write(f"{accession}\t{link}\n")
        logging.info(f"Download links saved for accession {accession}.")
    except Exception as e:
        logging.error(f"Error saving links for accession {accession}: {e}")

def main():
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Extract download links (AWS S3, EBI FTP, HTTP) for SRA/ENA accessions."
    )
    parser.add_argument("--accession", type=str, help="Single SRA/ENA accession.")
    parser.add_argument("--list", type=str, help="File containing a list of accessions.")
    args = parser.parse_args()

    accessions = []
    if args.accession:
        accessions.append(args.accession.strip())
    if args.list:
        if not os.path.exists(args.list):
            logging.error(f"List file {args.list} not found.")
            sys.exit(1)
        with open(args.list, "r") as f:
            for line in f:
                accession = line.strip()
                if accession:
                    accessions.append(accession)

    if not accessions:
        logging.error("No accession provided. Use --accession or --list.")
        sys.exit(1)

    for accession in accessions:
        logging.info(f"Processing accession: {accession}")
        links = get_sra_links(accession)
        if links:
            save_links(accession, links)
        else:
            logging.warning(f"Failed to retrieve links for accession {accession}")

    logging.info("Processing complete.")

if __name__ == "__main__":
    main()
