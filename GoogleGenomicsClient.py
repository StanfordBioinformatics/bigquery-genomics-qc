import argparse
import httplib2
import pprint
from apiclient.discovery import build
from collections import Counter
from oauth2client import tools
from oauth2client.client import flow_from_clientsecrets
from oauth2client.file import Storage
from oauth2client.tools import run_flow
import logging
from argparse import Namespace

class GoogleGenomicsClient(object):
    def __init__(self, client_secrets, project_number, dataset):
        self.client_secrets_path = client_secrets
        self.project_number = project_number
        self.service = self.setup(self.client_secrets_path)
        self.dataset = dataset

    def setup(self, client_secrets):
        storage = Storage('credentials.dat')
        credentials = storage.get()

        if credentials is None or credentials.invalid:
            flow = flow_from_clientsecrets(
                client_secrets,
                scope='https://www.googleapis.com/auth/genomics',
                message='You need to copy a client_secrets.json file into this directory, '
                'or pass in the --client_secrets_filename option to specify where '
                'one exists. See the README for more help.')

            # There's probably a better way to generate the 'flags' object.  Doing it this way for now.
            #parser2 = argparse.ArgumentParser(description=__doc__,
            #    formatter_class=argparse.RawDescriptionHelpFormatter,
            #    parents=[tools.argparser])
            #print parser2
            #parse2 = parser2.parse_args()
            #print parse2
            #parser2.add_argument('--client_secrets_filename',
            #    default=client_secrets,
            #    help='The filename of a client_secrets.json file from a '
            #         'Google "Client ID for native application" that '
            #         'has the Genomics API enabled.')

            flags = Namespace(client_secrets_filename=client_secrets, logging_level='INFO', auth_host_name='localhost',
                              auth_host_port=[8080, 8090], noauth_local_webserver=False)

            #flags = parser2.parse_args()

            credentials = run_flow(flow, storage, flags)
            # Create a genomics API service
        http = httplib2.Http()
        http = credentials.authorize(http)
        service = build('genomics', 'v1beta2', http=http)
        return service

    def list_datasets(self):
        request = self.service.datasets().list(projectNumber=self.project_number)
        response = self.execute(request)
        return response

    def get_call_set_id(self, call_set_name):
        request = self.service.callsets().search(
            body={'variantSetIds': [self.dataset], 'name': call_set_name}
        )
        response = self.execute(request)
        if "callSets" in response:
            for field in response["callSets"]:
                return field["id"]
        return None

    def delete_call_set(self, call_set_id):
        logging.debug("Deleting sample %s" % call_set_id)
        # todo: not executing until testing is complete
        request = self.service.callsets().delete(
            callSetId=call_set_id
        )
        response = self.execute(request)
        print response

    def get_variant_id(self, reference_name, start, end):
        request = self.service.variants().search(
            body={'variantSetIds': [self.dataset],
                  'referenceName': reference_name,
                  'start': start,
                  'end': end}
        )
        response = self.execute(request)
        if "variants" in response:
            for field in response["variants"]:
                return field["id"]
        return None

    def delete_variant(self, variant_id):
        # todo: not executing until testing is complete
        print "Deleting variant %s" % variant_id
        return None # todo turn turn this off
        #request = self.service.variants().delete(
        #    variantId=variant_id
        #)
        #response = self.execute(request)
        #print response

    def execute(self, request, retry=0):
        try:
            response = request.execute()
            if response is None and retry > 0:
                print "Request failed.  Retrying. %s" % request
                self.execute(request, retry-1)
            return response
        except Exception as e:
            logging.error("Request failed: %s" % e)


    def coverage_buckets(self, read_group, reference_name=None, start=None, end=None, ):
        request = self.service.readgroupsets().coveragebuckets().list(
            readGroupSetId=read_group,
            range_referenceName=reference_name,
            range_start=start,
            range_end=end,
        )
        response = self.execute(request, retry=2)
        #coverages = []
        if response is None:
            return None
        if "coverageBuckets" in response:
            return response["coverageBuckets"]
        return None

    def get_read_groups(self, dataset=None, sample_id=None):
        if dataset is None:
            dataset = self.dataset
        request = self.service.readgroupsets().search(
            body={'datasetIds': [dataset]}
        )
        response = self.execute(request)
        read_groups = []
        if "readGroupSets" in response:
            for field in response["readGroupSets"]:
                read_groups.append(field["id"])
        return read_groups