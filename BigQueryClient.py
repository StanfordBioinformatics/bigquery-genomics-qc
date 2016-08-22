import httplib2
from apiclient.discovery import build
from oauth2client.client import flow_from_clientsecrets
from oauth2client.file import Storage
from oauth2client import tools
import logging
from time import sleep
import json

class BigQueryClient(object):
    def __init__(self, client_secrets):
        self.client_secret_path = client_secrets

    def bigquery_setup(self):
        FLOW = flow_from_clientsecrets(self.client_secret_path,
                                   scope='https://www.googleapis.com/auth/bigquery')

        storage = Storage('bigquery_credentials.dat')
        credentials = storage.get()

        if credentials is None or credentials.invalid:
            credentials = tools.run_flow(FLOW, storage, tools.argparser.parse_args([]))

        http = httplib2.Http()
        http = credentials.authorize(http)

        bigquery_service = build('bigquery', 'v2', http=http)
        return bigquery_service

class BigQuery(object):
    def __init__(self, project_number, client_secrets, qc_dataset=None, project_name=None):
        client = BigQueryClient(client_secrets)
        self.service = client.bigquery_setup()
        self.project_number = project_number
        self.project_name = project_name
        self.qc_dataset = qc_dataset

    def run(self, query, query_name, table_output=False):
        logging.debug("Querying BigQuery")
        logging.debug(query)
        response = self.execute_query(query, query_name, table_output)
        if table_output is False:
            if self.check_response(response) is False:
                return False
        else:
            if self.check_insert(response) is True:
                job_id = self.get_job_id(response)
                return job_id
            else:
                return False
        result = self.parse_bq_response(response)
        return result

    # Run a query on BigQuery
    def execute_query(self, query, query_name, table_output=False):
        query_dict = {}
        if table_output is False:
            query_dict = {'query': query, 'timeoutMs': 1000000}
        else:
            query_dict = {
                'configuration': {
                    'query': {'query': query,
                        'allowLargeResults': 'true',
                        'destinationTable': {
                            'projectId': self.project_name,
                            'datasetId': self.qc_dataset,
                            'tableId': query_name
                        }
                    }
                }
            }

        query_request = self.service.jobs()
        try:
            query_response = None
            if table_output is False:
                query_response = query_request.query(projectId=self.project_number, body=query_dict).execute()
            else:
                query_response = query_request.insert(projectId=self.project_number, body=query_dict).execute()
            return query_response
        except Exception as e:
            logging.error("Unable to execute query: %s" % e)

    # Parse BigQuery response into list of dictionaries
    def parse_bq_response(self, response):
        fields = self.get_fields(response)
        result = []
        if 'rows' in response:
            for row in response['rows']:
                result_row = {}
                position = 0
                for field in row['f']:
                    result_row[fields[position]] = field['v']
                    position += 1
                result.append(result_row)
            return result
        return None

    # Get field names from BigQuery response
    def get_fields(self, response):
        fields = {}
        if 'schema' in response:
            position = 0
            for field in response['schema']['fields']:
                fields[position] = field['name']
                position += 1
            return fields
        return None

    # Check for valid response
    # todo add more checks
    def check_response(self, response):
        if response is None:
            logging.error("No Response)")
        if 'jobComplete' in response:
            if response['jobComplete'] is True:
                return True
            else:
                logging.error("Query timed out")
                exit(0)
        return False

    def get_job_id(self, response):
        job_reference =  response[u'jobReference']
        id = job_reference[u'jobId']
        return id

    def check_insert(self, response):
        status = response[u'status']
        if u'errors' in status:
            logging.error("BigQuery Error: %s" % status[u'errors'])
            return False
        return True

    def poll_job(self, job_id):
        state = 'RUNNING'
        while state == 'RUNNING':
            state = self.get_job_status(job_id)
            if state == 'RUNNING':
                sleep(5)
        return state

    def get_job_status(self, job_id):
        service = self.service.jobs()
        response = service.get(projectId=self.project_number, jobId=job_id).execute()
        if 'status' not in response:
            return False
        state = response['status']['state']
        return state