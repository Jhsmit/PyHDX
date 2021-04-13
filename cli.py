import argparse
from ipaddress import ip_address
from pyhdx.panel import serve
from pyhdx.panel.configurations import ConfigurationSettings

def main():
    parser = argparse.ArgumentParser(prog ='pyhdx',
                                     description ='PyHDX Launcher')
 
    parser.add_argument('serve', help= 'Runs PyHDX Dashboard')

    parser.add_argument('--cluster', help= 'Run with local cluster <ip> <port>', dest='cluster', nargs=2, metavar=('IP','PORT'))
    parser.add_argument('--cluster-auto', help='Runs dask locally', action='store_true',dest='cluster_auto')
  
    args = parser.parse_args()

    if (args.cluster):
        try:
            if (ip_address(args.cluster[0])) and int(args.cluster[1])>0 and int(args.cluster[1])<2**16: 
                ConfigurationSettings().update_cluster(args.cluster[0], args.cluster[1])
                print("Welcome to PyHDX server!")
                serve.run_main()                
        except ValueError:
            print("invalid IP address or Port given")
    elif args.cluster_auto:
            #yet to complete this
            print("will run cluster auto")            
    elif args.serve:
        print("Welcome to PyHDX server!")
        serve.run_main()
            

if __name__ == '__main__':
        main()
