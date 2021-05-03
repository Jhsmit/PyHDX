import argparse
from ipaddress import ip_address
from pyhdx.panel import serve
from pyhdx.panel.configurations import ConfigurationSettings
from dask.distributed import Client, LocalCluster
from pyhdx.support import verify_cluster


def main():
    parser = argparse.ArgumentParser(prog='pyhdx',
                                     description='PyHDX Launcher')
 
    parser.add_argument('serve', help="Runs PyHDX Dashboard")

    parser.add_argument('--cluster', help="Run with local cluster <ip> <port>", dest='cluster', nargs=2,
                        metavar=('IP', 'PORT'))
    args = parser.parse_args()

    if args.cluster:
        try:
            if (ip_address(args.cluster[0])) and 0 < int(args.cluster[1]) < 2 ** 16:
                ConfigurationSettings().update_cluster(args.cluster[0], args.cluster[1])
        except ValueError:
            print("Invalid IP address or Port given")
            return
        cluster = ':'.join(args.cluster)
        if not verify_cluster(cluster):
            return

    else:
        cluster = ConfigurationSettings().load_cluster()
        if not verify_cluster(cluster):
            # Start a new local cluster if none is specified
            local_cluster = LocalCluster()  #todo cluster config in configuration file
            _, ip, port = local_cluster.scheduler_address.split(':')
            ip = ip.strip('/')
            ConfigurationSettings().update_cluster(ip, port)

    if args.serve:
        serve.run_main()
            

if __name__ == '__main__':
    main()
