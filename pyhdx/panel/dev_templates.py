from jinja2 import Environment, FileSystemLoader
import os

pth = os.path.dirname(__file__)

env = Environment(loader=FileSystemLoader(pth))



template = env.get_template('golden.html')