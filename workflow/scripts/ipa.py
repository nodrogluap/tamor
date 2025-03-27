# This script defines functions for conducting Ingenuity Pathway Analysis (IPA) on a dataset
# and checking the status and results of the analysis.

# Import necessary libraries
from requests import get
import http.server
from urllib.parse import urlparse, parse_qs
from requests_oauthlib import OAuth2Session
import secrets
import hashlib
import base64
import webbrowser
import threading
import pandas as pd
import requests
import time
import subprocess

# Install dependencies


def install_dependencies():
    dependencies = [
        'requests',
        'http.server',
        'requests_oauthlib',
        'pandas'
    ]

    for dependency in dependencies:
        subprocess.call(['pip3', 'install', dependency])

# Define a function to start the OAuth callback server in a separate thread


def ipa_init():
    server_thread = threading.Thread(target=run_server)
    server_thread.daemon = True
    server_thread.start()

# Define a custom HTTP request handler for handling OAuth callbacks


class OAuthCallbackHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        # Parse the URL and query parameters
        parsed_url = urlparse(self.path)
        query_params = parse_qs(parsed_url.query)

        # Check if the request is for OAuth authorization
        if parsed_url.path == '/' and 'code' in query_params:
            # Extract authorization code and state from the query parameters
            auth_code = query_params['code'][0]
            state_code = query_params['state'][0]

            # Respond to the request
            self.send_response(200)
            self.send_header('Content-type', 'text/html')
            self.end_headers()
            self.wfile.write(
                b"Authorization received. You can close this window now.")

            # Pass the authorization code and state back to the main thread
            OAuthCallbackHandler.authorization_code = auth_code
            OAuthCallbackHandler.state_code = state_code
        else:
            # Respond with a 404 error for unknown paths or missing authorization code
            self.send_response(404)
            self.end_headers()
            self.wfile.write(b"Not Found")

# Function to run the HTTP server for handling OAuth callbacks


def run_server():
    server_address = ('', 8000)
    httpd = http.server.HTTPServer(server_address, OAuthCallbackHandler)
    #print('Server started at http://localhost:8000')
    httpd.serve_forever()

# Function to perform IPA analysis


def ipa_analyze(ipa, dataset_file, projectname, geneidtype, obs_names, measurement_types, cutoffs, datasetname=None, referenceset="dataset"):
    # Load dataset from file
    dataset = pd.read_csv(dataset_file, sep="\t", header=None)

    # Determine the dimensions of the dataset
    num_obs = len(obs_names)
    num_measurements = len(measurement_types)

    # Check if the dataset structure matches expectations
    if dataset.shape[1] != 1 + num_obs * num_measurements:
        print("Unexpected dataset structure")
        return

    # Derive dataset name from file name if not specified
    if datasetname is None:
        tokens = dataset_file.split("/")
        last_token = tokens[-1]
        datasetname = last_token.split(".")[0]

    post_data_list = []

    # Construct post data for IPA analysis
    post_data = "&".join([
        f"applicationname={ipa['application_name']}",
        f"projectname={projectname}",
        "ipaview=none",
        f"datasetname={datasetname}",
        f"analysisname={datasetname}",
        f"referenceset={referenceset}",
        f"geneidtype={geneidtype}",
        f"genecolname={dataset.iloc[0, 0]}"
    ])
    post_data_list.append(post_data)

    # Append observation names to post data
    for i, obs_name in enumerate(obs_names):
        post_data_list.append(f"&obs{i+1}name={obs_name}")

    # Append measurement types to post data
    for i, measurement_type in enumerate(measurement_types):
        if i == 0:
            post_data_list.append(f"&expvaltype={measurement_type}")
        else:
            post_data_list.append(f"&expvaltype{i+1}={measurement_type}")

    # Append measurement column headers to post data
    for i, obs_name in enumerate(obs_names):
        for j, measurement_type in enumerate(measurement_types):
            dataset_col = (i * num_measurements) + j + 1
            if j == 0:
                post_data_list.append(
                    f"&{'obs'+str(i+1) if i > 0 else ''}expvalname={dataset.iloc[0, dataset_col]}")
            else:
                post_data_list.append(
                    f"&{'obs'+str(i+1) if i > 0 else ''}expval{j+1}name={dataset.iloc[0, dataset_col]}")

    # Append cutoffs to post data
    for i, cutoff in enumerate(cutoffs):
        if not pd.isna(cutoff):
            if i == 0:
                post_data_list.append(f"&cutoff={cutoff}")
            else:
                post_data_list.append(f"&cutoff{i+1}={cutoff}")

    # Append genes and measurements to post data
    for i in range(1, dataset.shape[0]):
        post_data = ""
        for j in range(dataset.shape[1]):
            if j == 0:
                param = "geneid"
            else:
                obs_index = (j - 1) // num_measurements + 1
                type_index = (j - 1) % num_measurements + 1
                if type_index == 1:
                    param = f"expvalue"
                else:
                    param = f"expval{type_index}"
            value = str(dataset.iloc[i, j]) if not pd.isna(
                dataset.iloc[i, j]) else "NaN"
            post_data += f"&{param}={value}"
        post_data_list.append(post_data)

    post_data = "".join(post_data_list)

    # Make POST request to perform IPA analysis
    url = f"https://{ipa['host']}/pa/api/v2/multiobsanalysis"
    headers = {'Authorization': f"Bearer {ipa['access_token']}",
               'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post(url, headers=headers, data=post_data)

    # Extract analysis IDs from response
    analysis_ids = response.text.split(",")
    return analysis_ids

# Function to request OAuth authorization code


def ipa_login(authorization_base_url, token_url, client_id='1571511054-1646124475-167497993-BQfZBb', scope=None, state=None, application_name='PythonAPI', redirect_uri='http://localhost:8000'):
    # Generate code verifier and code challenge
    code_verifier = secrets.token_urlsafe(100)
    code_challenge = base64.urlsafe_b64encode(hashlib.sha256(
        code_verifier.encode()).digest()).decode().rstrip('=')

    # Create OAuth session without code_verifier
    oauth = OAuth2Session(
        client_id, redirect_uri=redirect_uri, scope=scope, state=state)

    # Construct authorization URL with code challenge and method
    authorization_url, _ = oauth.authorization_url(
        authorization_base_url,
        code_challenge=code_challenge,
        code_challenge_method='S256'  # Specify the code challenge method
    )

    # Open authorization URL in browser for user to authorize
    webbrowser.open(authorization_url)

    # Wait for authorization response
    while not hasattr(OAuthCallbackHandler, 'authorization_code'):
        pass

    # Use the authorization code and state from the callback
    auth_code = OAuthCallbackHandler.authorization_code
    state_code = OAuthCallbackHandler.state_code

    # Exchange authorization code for access token
    token = oauth.fetch_token(
        token_url, code=auth_code, code_verifier=code_verifier, state=state_code)

    # Extract access token and other required information
    access_token = token['access_token']
    cookie_file = token.get('cookieFile', None)

    # Define the dictionary to return
    oauth_info = {
        'host': 'analysis.ingenuity.com',
        'cookieFile': cookie_file,
        'application_name': application_name,
        'access_token': access_token
    }
    return oauth_info


def ipa_check_status(ipa, analysis_ids, ping_interval=30, max_attempts=100):
    # Check if running in a Jupyter notebook
    try:
        from IPython import get_ipython
        if 'IPKernelApp' not in get_ipython().config:
            # If not in a notebook, disable notebook-specific features
            from tqdm import tqdm_notebook as tqdm
        else:
            # If in a notebook, use tqdm
            from tqdm import tqdm
    except:
        # If unable to import IPython, assume running in a script
        from tqdm import tqdm

    if not isinstance(analysis_ids, list):
        analysis_ids = [analysis_ids]  # Ensure analysis_ids is a list

    num_analyses = len(analysis_ids)
    analyze_result = f"Submitted {num_analyses} {'analysis' if num_analyses == 1 else 'analyses'}: {', '.join(map(str, analysis_ids))}"
    print(analyze_result)

    attempt = 0
    analysis_status = []

    while attempt < max_attempts:
        print(f"Status check loop (attempt {attempt+1} of {max_attempts})")

        for analysis_id in tqdm(analysis_ids, desc="Checking status"):
            status = check_status(ipa, analysis_id)

            # Interpret the status
            if status in ['3', '4', '5']:  # Convert to string for comparison
                print(
                    f"DONE: Analysis {analysis_id} finished with status: {status}")
                analysis_status.append((analysis_id, status))
            else:
                # Analysis is not done yet
                print(f"Analysis {analysis_id} still in progress")

        if len(analysis_status) == num_analyses:
            # All analyses have finished
            break
        else:
            attempt += 1
            if attempt < max_attempts:
                time.sleep(ping_interval)

    return analysis_status


def check_status(ipa, analysis_id):
    url = f"https://{ipa['host']}/pa/api/v2/analysisstatus?applicationname={ipa['application_name']}&analysisuid={analysis_id}"
    headers = {"Authorization": f"Bearer {ipa['access_token']}"}

    response = requests.get(url, headers=headers)

    status = response.text
    return status

# Function to retrieve analysis results from IPA


def ipa_results(analysis_ids, status, ipa):
    canonical_pathways_results = []
    upstream_regulators_results = []
    biofunctions_results = []

    # Loop through analysis IDs
    for analysis_id in analysis_ids:
        found = False
        for analysis_tuple in status:
            if analysis_tuple[0] == analysis_id:
                found = True
                status_code = analysis_tuple[1]
                if status_code == '3':
                    result = 'Succeeded'
                elif status_code == '4':
                    result = 'Failed'
                else:
                    result = 'Canceled'
                analysis_result = {
                    'Analysis ID': analysis_id, 'Status': result}

                # Retrieve top analysis scores for different entities
                if status_code == '3':
                    cp_scores = ipa_get_top_analysis_scores(ipa, analysis_id, ipa_analysis_entity()[
                                                            "CANONICAL_PATHWAY"])
                    if cp_scores:
                        for i, score in enumerate(cp_scores, start=1):
                            cp_result = {
                                'Analysis ID': analysis_id, 'Status': result}
                            cp_result.update(score)
                            canonical_pathways_results.append(cp_result)
                    else:
                        canonical_pathways_results.append({'Analysis ID': analysis_id, 'Status': result, 'pathways': 'Failed to retrieve scores', 'Pvalue': 'Failed to retrieve scores',
                                                          'Zscore': 'Failed to retrieve scores', 'OverlappingGenes': 'Failed to retrieve scores', 'OverlappingGenesCount': 'Failed to retrieve scores'})

                    ur_scores = ipa_get_top_analysis_scores(ipa, analysis_id, ipa_analysis_entity()[
                                                            "UPSTREAM_REGULATOR"])
                    if ur_scores:
                        for i, score in enumerate(ur_scores, start=1):
                            ur_result = {
                                'Analysis ID': analysis_id, 'Status': result}
                            ur_result.update(score)
                            upstream_regulators_results.append(ur_result)
                    else:
                        upstream_regulators_results.append({'Analysis ID': analysis_id, 'Status': result, 'entity': 'Failed to retrieve scores', 'Pvalue': 'Failed to retrieve scores',
                                                           'Zscore': 'Failed to retrieve scores', 'OverlappingGenes': 'Failed to retrieve scores', 'OverlappingGenesCount': 'Failed to retrieve scores'})

                    bf_scores = ipa_get_top_analysis_scores(ipa, analysis_id, ipa_analysis_entity()[
                                                            "BIOFUNCTION"])
                    if bf_scores:
                        for i, score in enumerate(bf_scores, start=1):
                            bf_result = {
                                'Analysis ID': analysis_id, 'Status': result}
                            bf_result.update(score)
                            biofunctions_results.append(bf_result)
                    else:
                        biofunctions_results.append({'Analysis ID': analysis_id, 'Status': result, 'entity': 'Failed to retrieve scores', 'Pvalue': 'Failed to retrieve scores',
                                                    'Zscore': 'Failed to retrieve scores', 'OverlappingGenes': 'Failed to retrieve scores', 'OverlappingGenesCount': 'Failed to retrieve scores'})

                break  # Stop searching for the analysis ID once found
        if not found:
            canonical_pathways_results.append({'Analysis ID': analysis_id, 'Status': f'Analysis {analysis_id} did not finish in time',
                                              'pathways': 'N/A', 'Pvalue': 'N/A', 'Zscore': 'N/A', 'OverlappingGenes': 'N/A', 'OverlappingGenesCount': 'N/A'})
            upstream_regulators_results.append({'Analysis ID': analysis_id, 'Status': f'Analysis {analysis_id} did not finish in time',
                                               'entity': 'N/A', 'Pvalue': 'N/A', 'Zscore': 'N/A', 'OverlappingGenes': 'N/A', 'OverlappingGenesCount': 'N/A'})
            biofunctions_results.append({'Analysis ID': analysis_id, 'Status': f'Analysis {analysis_id} did not finish in time',
                                        'entity': 'N/A', 'Pvalue': 'N/A', 'Zscore': 'N/A', 'OverlappingGenes': 'N/A', 'OverlappingGenesCount': 'N/A'})

    return pd.DataFrame(canonical_pathways_results), pd.DataFrame(upstream_regulators_results), pd.DataFrame(biofunctions_results)

# Function to define IPA analysis entities


def ipa_analysis_entity():
    return {
        "CANONICAL_PATHWAY": "CANONICAL_PATHWAY",
        "UPSTREAM_REGULATOR": "UPSTREAM_REGULATOR",
        "BIOFUNCTION": "BIOFUNCTION"
    }


def ipa_get_analysis_result(ipa, analysis_id, entity):
    result_url = {
        ipa_analysis_entity()["CANONICAL_PATHWAY"]: "allCanonicalPathways",
        ipa_analysis_entity()["UPSTREAM_REGULATOR"]: "allUpstreamRegulators",
        ipa_analysis_entity()["BIOFUNCTION"]: "allBioFunctions"
    }.get(entity)

    if not result_url:
        print("Did not recognize analysis result type:", entity)
        return None

    url = f"https://{ipa['host']}/pa/ipa/analysisResults/{result_url}/{ipa['application_name']}/{analysis_id}"
    resp = requests.get(
        url, headers={"Authorization": f"Bearer {ipa['access_token']}"})

    if resp.status_code == 200:
        return resp.json()
    else:
        print("Failed to retrieve analysis result")
        return None


def ipa_get_top_analysis_scores(ipa, analysis_id, entity):
    json_data = ipa_get_analysis_result(ipa, analysis_id, entity)
    headers = {"CANONICAL_PATHWAY": "pathways",
               "UPSTREAM_REGULATOR": "regulators",
               "BIOFUNCTION": "functions"}

    top_scores = []

    for pathway_data in json_data.get(headers.get(entity), []):
        pathway_info = {}

        # Dynamically learn keys from the JSON data
        keys = pathway_data.keys()

        for key in keys:
            value = pathway_data[key]
            if isinstance(value, str):
                # Clean up name field if necessary
                if key == 'name':
                    value = value.replace("+", " ")
            elif isinstance(value, (int, float)):
                # Convert numerical values to float
                value = float(value)

            pathway_info[key.capitalize()] = value

        top_scores.append(pathway_info)

    # Sort by Pvalue with NaN at the end
        top_scores.sort(key=lambda x: (float('inf') if x.get("Pvalue") is None or x.get(
            "Pvalue") == "NaN" else float(x.get("Pvalue"))), reverse=False)

    return top_scores

