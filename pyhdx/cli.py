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
    parser.add_argument('--scheduler_address', help="Run with local cluster <ip>:<port>")
    args = parser.parse_args()

    cfg = ConfigurationSettings()

    if args.scheduler_address:
        ip, port = args.scheduler_address.split(':')
        if not ip_address(ip):
            print('Invalid IP Address')
            return
        elif not 0 <= int(port) < 2**16:
            print('Invalid port, must be 0-65535')
            return
        cfg.set('cluster', 'scheduler_address', args.scheduler_address)

    scheduler_address = cfg.get('cluster', 'scheduler_address')
    if not verify_cluster(scheduler_address):
        # Start a new local cluster if none is found
        client = default_cluster()
        _, ip, port = client.scheduler_address.split(':')
        ip = ip.strip('/')
        scheduler_address = f"{ip}:{port}"
        print(f"Started new Dask LocalCluster at {scheduler_address}")

    if args.serve:
        serve.run_main()
            

if __name__ == '__main__':
    import sys
    sys.argv.append('serve')
    sys.argv.append('--scheduler_address')
    sys.argv.append('127.0.0.1:53270')

    main()
