import configparser
from pathlib import Path

default_ip = '127.0.0.1'
default_port = '52123'

directory = Path(__file__).parent


class ConfigurationSettings(object):
    def load_config(self):
        """
        This method reads the configuration file
        """
        config = configparser.ConfigParser()
        config.read(self.config_file)

        return config

    def __init__(self, config_file='config.ini'):
        """
        Parameters
        ----------
        config_file : pass any *.ini file as argument to load and update, by default reads config.ini in pyhdx/panel
        """
        self.config_file = directory / config_file
        self.config = self.load_config()

    def load_cluster(self):
        """
        This method will load the dask server host IP and Port values from the configuration file
        Returns
        -------
        <ip>:<port> : :obj:`str`
            Default value - '127.0.0.1:52123'
        """
        if not self.config.has_section('cluster'):
            self.update_cluster(default_ip, default_port)
        return str(self.config.get('cluster','ip'))+":"+str(self.config.get('cluster','port'))

    def update_cluster(self, ip, port):
        """
        This method is to update the dask server host IP and Port values in the configuration file.
        Parameters
        ----------
        ip : :obj:`str` valid IPv4 address
        port : :obj:`str` valid port address less than 65535
        """
        if not self.config.has_section('cluster'):
            self.config.add_section('cluster')
        self.config.set('cluster', 'ip', ip)
        self.config.set('cluster', 'port', port)
        self.update_config()

    def update_config(self):
        """
        This method is to updates the configuration file.
        """
        with open(self.config_file, 'w') as config_file:
            self.config.write(config_file)
