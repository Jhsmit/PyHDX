import argparse
from ipaddress import ip_address
from pyhdx.web import serve
from pyhdx.config import ConfigurationSettings
from pyhdx.local_cluster import verify_cluster, default_cluster


# todo add check to see if the web module requirements are installed

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
        cfg.cluster = ':'.join(args.cluster)

    cluster = cfg.cluster
    if not verify_cluster(cluster):
        # Start a new local cluster if none is found
        client = default_cluster()
        _, ip, port = client.scheduler_address.split(':')
        ip = ip.strip('/')
        cluster = f"{ip}:{port}"
        print(f"Started new Dask LocalCluster at {cluster}")
        cfg.cluster = cluster

    if args.serve:
        serve.run_main()
            

if __name__ == '__main__':
    # import sys
    # sys.argv.append('serve')
    # sys.argv.append('--cluster')
    # sys.argv.append('127.0.0.1')
    # sys.argv.append('53270')
    main()
