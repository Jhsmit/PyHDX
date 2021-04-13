import os
import configparser
from pathlib import Path

default_ip = '127.0.0.1'
default_port = '52123'

directory = Path(__file__).parent

class ConfigurationSettings(object):
    def load_config(self):
        self.config = configparser.ConfigParser()
        self.config.read(self.config_file)

    def __init__(self, config_file='config.ini'):
        self.config_file = directory / config_file
        self.load_config()

    def load_cluster(self):
        if not self.config.has_section('cluster'):
            self.update_cluster(default_ip, default_port)
        return str(self.config.get('cluster','ip'))+":"+str(self.config.get('cluster','port'))

    def update_cluster(self, ip, port):
        if not self.config.has_section('cluster'):
            self.config.add_section('cluster')
        self.config.set('cluster', 'ip', ip)
        self.config.set('cluster', 'port', port)
        self.update_config()

    def update_config(self):
        with open(self.config_file, 'w') as config_file:
            self.config.write(config_file)
