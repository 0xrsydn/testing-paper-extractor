from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
import time
import os
import concurrent.futures
import pandas as pd
import autoit
from dotenv import load_dotenv

def main():
    
    configs = get_env()
    
    code = preprocess(configs)
    
    if code == 1:
        print("Please copy the files to another location and clear the download directory and try again.")
        return
    
    tokens = [120000, 96000, 64000, 32000]
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = [executor.submit(exec, configs, token) for token in tokens]
        for f in concurrent.futures.as_completed(results):
            print(f.result())
    
    rename_files(configs)
    merge_files(configs)
    upload_to_one_drive(configs)
    
    print('All Done!')
    
def get_env():
    load_dotenv()
    return {
        'paper_name': os.getenv('PAPER_NAME'),
        'hf_url': os.getenv('HF_URL'),
        'url': os.getenv('URL'),
        'extraction_mode': os.getenv('EXTRACTION_MODE'),
        'base_dir': os.getenv('BASE_DIR'),
        'one_drive_url': os.getenv('ONE_DRIVE_URL'),
    }
    
def preprocess(conf):
    print("Preprocessing...")
    print("Removing old files...")
    
    flag = False
    
    
    val = check_if_base_dir_exists(conf)
    
    if val == 0:
        return 2
    
    check_if_chrome_driver_exists()
    
    base_dir = conf['base_dir']
    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)
        if os.path.isdir(subdir_path):
            remove_files(subdir_path, flag)
            
    remove_files(base_dir, flag)
    
    return 0

def check_if_chrome_driver_exists():
    if not os.path.exists('chromedriver.exe'):
        print("Please download the Chrome WebDriver from https://sites.google.com/a/chromium.org/chromedriver/downloads and place it in the same directory as this script.")
        return 2
    return 0

def check_if_base_dir_exists(configs):
    base_dir = configs['base_dir']
    if not os.path.exists(base_dir):
        print("Creating directories...")
        os.makedirs(base_dir)
        for token in ['120000', '96000', '64000', '32000']:
            os.makedirs(os.path.join(base_dir, token))
        return 0
    return 1

def remove_files(directory, flag):
    files = os.listdir(directory)
    for file in files:
        if file.endswith('.xlsx'):
            while not flag:
                ans = input("Old files found. Do you want to remove them? (y/n): ")
                if ans.lower() == 'y':
                    flag = True
                elif ans.lower() == 'n':
                    return 1
                else:
                    print("Invalid input. Please enter 'y' or 'n'.")
            os.remove(os.path.join(directory, file))
            print(f"Removed {file}")
    return 0
    
    
def exec(configs, token):
    options = webdriver.ChromeOptions()
    prefs = {
        "download.default_directory": os.path.abspath(f"./downloads/{token}"),
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True
    }
    options.add_experimental_option("prefs", prefs)
    # options.add_argument("--headless")
    driver = webdriver.Chrome(options=options)
    driver.get(configs['url'])
    
    driver.implicitly_wait(5)
    
    try:
        WebDriverWait(driver, 100000).until(EC.text_to_be_present_in_element((By.ID, 'nutrigenme-paper-extraction'), "NutriGenMe - Paper Extraction"))
    except TimeoutException:
        print("Text NutriGenMe - Paper Extraction is not present on the page within the specified time.")
    
    abs_paper_path = os.path.abspath(configs['paper_name'])
    buttons = driver.find_elements(By.XPATH, '//*[@data-testid="stFileUploaderDropzoneInput"]')
    buttons[0].send_keys(abs_paper_path)
    
    driver.implicitly_wait(5)
    
    toggle_checkbox = driver.find_elements(By.XPATH, '//*[@data-baseweb="checkbox"]')
    if configs['extraction_mode'] == 'text':
        toggle_checkbox[1].click()
    else:
        toggle_checkbox[0].click()
        
    driver.implicitly_wait(5)
    
    dropdown = driver.find_element(By.XPATH, '//*[@data-baseweb="select"]')
    dropdown.click()
    
    driver.implicitly_wait(2)
    
    if(token == 120000):
        token_id = "0"
    elif(token == 96000):
        token_id = "1"
    elif(token == 64000):
        token_id = "2"
    else:
        token_id = "3"
    
    selected_option = driver.find_element(By.ID, f'bui6val-{token_id}')
    selected_option.click()
    
    start_button = driver.find_elements(By.XPATH, '//*[@data-testid="baseButton-secondary"]')
    start_button[1].click()
    
    expected_text = "Gene and SNPs succesfully collected."
    
    try:
        WebDriverWait(driver, 100000).until(EC.text_to_be_present_in_element((By.XPATH, "//*"), expected_text))
    except TimeoutException:
        print("Text '{}' is not present on the page within the specified time.".format(expected_text))
    
    try:
        WebDriverWait(driver, 120).until(EC.text_to_be_present_in_element((By.XPATH, "//*"), "Save Result"))
    except TimeoutException:
        print("Text '{}' is not present on the page within the specified time.".format("Save Result"))
    
    get_values = driver.find_elements(By.TAG_NAME, 'code')
    print('Entities Extraction Done', get_values[0].text)
    print('Generating Summary Done', get_values[1].text)
    print('Table Extraction Done', get_values[2].text)
    
    save_button = driver.find_elements(By.XPATH, '//*[@data-testid="baseButton-secondary"]')
    save_button[2].click()
    
    print("Please wait...")
    print(f"Downloading the file with token {token}...")
    time.sleep(10)
    
    driver.close()
    
def rename_files(conf):
    base_dir = conf['base_dir']
    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)
        if os.path.isdir(subdir_path):
            files = os.listdir(subdir_path)
            if len(files) == 1 and files[0].endswith('.xlsx'):
                old_file_path = os.path.join(subdir_path, files[0])
                new_file_path = os.path.join(subdir_path, f"{subdir}.xlsx")
                os.rename(old_file_path, new_file_path)
                print(f"Renamed {files[0]} to {subdir}.xlsx")
                
def merge_files(conf):
    base_dir = conf['base_dir']
    merged_file_path = os.path.join(base_dir, conf['paper_name'].replace('.pdf', '.xlsx'))
    with pd.ExcelWriter(merged_file_path) as writer:
        sorted_subdirs = sorted(os.listdir(base_dir), key=lambda x: int(x) if x.isdigit() else float('inf'))
        for subdir in sorted_subdirs:
            subdir_path = os.path.join(base_dir, subdir)
            if os.path.isdir(subdir_path):
                files = os.listdir(subdir_path)
                if len(files) == 1 and files[0].endswith('.xlsx'):
                    df = pd.read_excel(os.path.join(subdir_path, files[0]))
                    df.to_excel(writer, sheet_name=subdir, index=False)
                
def upload_to_one_drive(configs):
    options = webdriver.ChromeOptions()
    # options.add_argument("--headless")
    driver = webdriver.Chrome(options=options)
    driver.get(configs['one_drive_url'])

    try:
        upload_button = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.XPATH, "//button[@title='Upload files from your computer to this location']")))
        upload_button.click()
        time.sleep(2)
        upload_folder_button = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.XPATH, "//button[@data-automationid='uploadFileCommand']")))
        upload_folder_button.click()
        time.sleep(2)
        autoit.control_send("Open", "Edit1", os.path.abspath(f"./downloads/{configs['paper_name'].replace('.pdf', '.xlsx')}"))
        autoit.control_send("Open", "Edit1", "{ENTER}")
        
        print("Please wait...")
        print("Uploading the file...")
        time.sleep(10)
        
        print("File uploaded successfully!")
    except TimeoutException:
        print("Failed to upload the file.")

    driver.quit()

if __name__ == '__main__':
    main()
