import h5py
import xml.etree.ElementTree as ET

from zipfile import ZipFile


def ods(name):
    with ZipFile(name) as z:
        with z.open('content.xml') as f:
            tree = ET.parse(f)
    root = tree.getroot()
