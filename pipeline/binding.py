"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""
from pipeline import PipelineElement
import urllib2, urllib
import pandas as pd
from StringIO import StringIO

class IEDBMHCBinding(PipelineElement):

  def __init__(self, alleles=[], name="IEDB-MHC-Binding", method='recommended', lengths = [8,9,10,11], url='http://tools.iedb.org/tools_api/mhci/'):
    self.name = name
    self._method = method
    self._lengths = lengths
    self._url = url
    self._alleles = ",".join(alleles)

  def query_iedb(self, sequence):
    request_values = {

        "method" : self._method,
        "length" : ",".join(str(l) for l in self._lengths),
        "sequence_text" : sequence,
        "allele" : self._alleles,
    }

    print "Calling iedb with", sequence, self._alleles
    try:
      data = urllib.urlencode(request_values)
      req = urllib2.Request(self._url, data)
      response = urllib2.urlopen(req).read()

      return pd.read_csv(StringIO(response), sep='\t', na_values=['-'])
    except:
      print "Connection error: Failed on sequence", sequence
      return pd.DataFrame()

  def apply(self,data):
    responses = []
    for peptide in data:
       responses += [self.query_iedb(peptide)]
    

    return pd.concat(responses)



if __name__ == '__main__':
  iedb = IEDBMHCBinding(alleles=['HLA-C*12:03', 'HLA-C*12:02'])
  print iedb.query_iedb("APHHSGVYPVNVQLYEAWKKV")
