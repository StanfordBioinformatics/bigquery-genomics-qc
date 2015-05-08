import sys
sys.path.append('/Users/gmcinnes/GitHub/mvp_aaa_codelabs/qc/qc_pipeline')
#import GoogleGenomicsClient
from GoogleGenomicsClient import GoogleGenomicsClient
from config import Config

client_secret_path = Config.CLIENT_SECRETS
project_number = Config.PROJECT_NUMBER
dataset = Config.DATA_SET
g = GoogleGenomicsClient(client_secrets=client_secret_path, project_number=project_number, dataset=dataset)


response = g.get_call_set_id(sample_name='LP6005038-DNA_B01')
print response
