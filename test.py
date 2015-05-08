import sys
sys.path.append('/Users/gmcinnes/GitHub/mvp_aaa_codelabs/qc/qc_pipeline')
#import GoogleGenomicsClient
from GoogleGenomicsClient import GoogleGenomicsClient
from config import Config

client_secret_path = Config.CLIENT_SECRETS
project_number = Config.PROJECT_NUMBER
dataset = Config.DATASET
g = GoogleGenomicsClient(client_secrets=client_secret_path, project_number=project_number, dataset=dataset)


#response = g.get_call_set_id(call_set_name='LP6005038-DNA_B01')
#print response

response = g.delete_variant(variant_id='CNimi86fy8umdhIEY2hyMRjQ-C0gvOCFm_HK44B-YZ')
print response
'''
{
 "start": "752565",
 "variantSetIds": [
  "8524520633659478872"
 ],
 "referenceName": "chr1",
 "maxCalls": 1,
 "end": "752566"
}
'''