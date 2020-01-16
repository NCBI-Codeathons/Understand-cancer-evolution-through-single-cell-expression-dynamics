 import os
import sys
import tarfile
import pandas as pd
import gzip
from io import BytesIO
from io import StringIO
import numpy as np
import csv
import gzip
import GEOparse
import re

#from datasets import geoREAD
class geoREAD:

    def __init__(self, gse_id, data_dir):
        #Data directory with raw tar files from GEO
        self.data_dir = self._mkdir(data_dir)
        #Initialize the GSE id file of the dataset
        self.gse_id = self._getData(gse_id)
        #GSE tar file name
        self.gse_file = os.path.join(self.data_dir,gse_id)
        #Gsm files
        self.gsm_files = self._getFiles(gse_id)

    def _mkdir(self, data_dir):
        if not os.path.exists(data_dir):
            ans = input('Do you want to create the dir %s (y|n):' % (data_dir))
            if ans.lower().strip()[0] == 'y':
                os.mkdir(data_dir)

            else:
                raise Exception("Make folder when ready...")
                return -1

        return data_dir

    def _getData(self, gse_id):
        """Add functionality for tar.gz files
        """
        gse_file = os.path.join(self.data_dir,gse_id)
        if not os.path.isfile(gse_file):
            raise Exception('Raw {} data missing. Download from GEO and put in {}'\
            .format(gse_id,self.data_dir))
            return -1

        else:
            return gse_id

    def _getFiles(self, gse_id):
        with tarfile.open(self.gse_file, 'r') as tar:
            gsm_files = tar.getnames()

        return gsm_files

    def readGSM(self, gsm_filename):
        with tarfile.open(self.gse_file, 'r') as tar:
            for member in [tar.getmember(gsm_filename)]:
                file = tar.extractfile(member)
                with gzip.open(file, 'rt', newline='') as f:
                    data = pd.read_csv(f)

        return data

#Debugging and testing
if __name__ == "__main__":

    #Initialize a data object
    geo = geoREAD('GSE118828_RAW.tar','./data_set/')
    scdata = geo.readGSM(geo.gsm_files[1])
    print(scdata.shape)
    print(scdata.iloc[:5,:5])
