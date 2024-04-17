# NutrigenMe Paper Extractor Automation Tool

## Description
This is a CLI tool is designed to extract the text from a PDF file and save it to a xlsx file automatically then upload it to the designated Microsoft OneDrive folder.

## Installation
1. Clone the repository
2. Make an environment using the following command:
```bash
python -m venv env
```
3. Activate the environment using the following command:
Unix
```bash
source env/bin/activate
```
Windows
```bash
env\Scripts\activate
```
4. Install the required packages using the following command:
```bash
pip install -r requirements.txt
```
5. Put the PDF file in the same directory as the script.
6. Make sure to have downloaded the chromedriver from [here](https://chromedriver.chromium.org/downloads) and put it in the same directory as the script.
7. Create and configure the variables in a .env file with the following format:
```env
PAPER_NAME=<YOUR_PDF_NAME_WITH_EXTENSION> -> e.g. NutriGenMePE-2.0.pdf
HF_URL=https://huggingface.co/spaces/KalbeDigitalLab/NutriGenMePE-2.0
URL=https://kalbedigitallab-nutrigenmepe-2-0.hf.space/
EXTRACTION_MODE=<EXTRACTION_MODE> -> e.g. text or table
BASE_DIR=downloads
ONE_DRIVE_URL=<ONE_DRIVE_URL_WITH_NO_AUTH_REQUIRED>
```
8.  Run the script using the following command:
```bash
python bot.py
```
9.  Wait for the script to finish and check the designated OneDrive folder for the xlsx file.