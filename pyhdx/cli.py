import argparse
from ipaddress import ip_address
from pyhdx.panel import serve
from pyhdx.panel.config import ConfigurationSettings
from dask.distributed import Client, LocalCluster
from pyhdx.support import verify_cluster


def main():
    parser = argparse.ArgumentParser(prog='pyhdx',
                                     description='PyHDX Launcher')
 
    parser.add_argument('serve', help="Runs PyHDX Dashboard")

    parser.add_argument('--cluster', help="Run with local cluster <ip> <port>", dest='cluster', nargs=2,
                        metavar=('IP', 'PORT'))
    args = parser.parse_args()

    cfg = ConfigurationSettings()

    if args.cluster:
        if not ip_address(args.cluster[0]):
            print('Invalid IP Address')
            return
        elif not 0 < int(args.cluster[1]) < 2**16:
            print('Invalid port, must be 0-65535')
            return
        elif not verify_cluster(':'.join(args.cluster)):
            print("No Dask scheduler found at given address")
            return
        cfg.cluster = ':'.join(args.cluster)

    else:
        cluster = ConfigurationSettings().cluster
        if not verify_cluster(cluster):
            # Start a new local cluster if none is specified
            local_cluster = LocalCluster()  #todo cluster config in configuration file
            _, ip, port = local_cluster.scheduler_address.split(':')
            ip = ip.strip('/')
            cluster = f"{ip}:{port}"
            print(f"Started new Dask LocalCluster at {cluster}")
            cfg.cluster = cluster

    if args.serve:
        serve.run_main()
            

if __name__ == '__main__':
    import sys
    sys.argv.append('serve')
    main()
