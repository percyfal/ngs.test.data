"""ngs test data module"""
import os
import yaml

__import__('pkg_resources').declare_namespace(__name__)

MODULEDIR=os.path.dirname(__file__)

# Configuration file
with open (os.path.join(os.path.dirname(__file__), os.pardir, "config", "config.yaml")) as fh:
    config = yaml.load(fh)
